#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robust RNA-seq Pipeline with Resume Capability
Created for overnight processing with error recovery
"""

import os
import json
import time
import logging
import pysam
import subprocess
import pandas as pd
from pathlib import Path
from datetime import datetime
from plotting_reads_from_sorted_BAM_file import plotting_sorted_BAM

# Setup logging
def setup_logging():
    """Setup comprehensive logging"""
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"pipeline_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # Also print to console
        ]
    )
    return logging.getLogger(__name__)

class RobustRNASeqPipeline:
    def __init__(self, config_file="pipeline_config.json"):
        self.logger = setup_logging()
        self.config = self.load_config(config_file)
        self.progress_file = "pipeline_progress.json"
        self.progress = self.load_progress()
        
    def load_config(self, config_file):
        """Load configuration with defaults"""
        default_config = {
            "tools": {
                "bowtie2_build": "bowtie2-build",
                "bowtie2": "bowtie2", 
                "fastq_dump": "fastq-dump"
            },
            "directories": {
                "output": "sorted_bam_files",
                "temp": "temp_files",
                "home": str(Path.home())
            },
            "parameters": {
                "bowtie2_threads": "4",
                "max_retries": 3,
                "retry_delay": 60,  # seconds
                "download_timeout": 1800  # 30 minutes
            }
        }
        
        if os.path.exists(config_file):
            with open(config_file, 'r') as f:
                user_config = json.load(f)
                # Merge with defaults
                for section in default_config:
                    if section in user_config:
                        default_config[section].update(user_config[section])
        else:
            # Create default config file
            with open(config_file, 'w') as f:
                json.dump(default_config, f, indent=2)
            self.logger.info(f"Created default config file: {config_file}")
            
        return default_config
    
    def load_progress(self):
        """Load progress from previous run"""
        if os.path.exists(self.progress_file):
            with open(self.progress_file, 'r') as f:
                progress = json.load(f)
                self.logger.info(f"Resuming from previous run. Completed: {len(progress.get('completed', []))} files")
                return progress
        return {"completed": [], "failed": [], "current_step": "init"}
    
    def save_progress(self):
        """Save current progress"""
        with open(self.progress_file, 'w') as f:
            json.dump(self.progress, f, indent=2)
    
    def run_command_with_retry(self, command, description, max_retries=None):
        """Run command with retry logic"""
        if max_retries is None:
            max_retries = self.config["parameters"]["max_retries"]
            
        for attempt in range(max_retries):
            try:
                self.logger.info(f"Attempt {attempt + 1}/{max_retries}: {description}")
                result = subprocess.run(
                    command, 
                    check=True, 
                    capture_output=True, 
                    text=True,
                    timeout=self.config["parameters"]["download_timeout"]
                )
                self.logger.info(f"Success: {description}")
                return result
                
            except subprocess.TimeoutExpired:
                self.logger.warning(f"Timeout for {description} (attempt {attempt + 1})")
                if attempt < max_retries - 1:
                    time.sleep(self.config["parameters"]["retry_delay"])
                    
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Command failed: {description}")
                self.logger.error(f"Error: {e.stderr}")
                if attempt < max_retries - 1:
                    time.sleep(self.config["parameters"]["retry_delay"])
                    
            except Exception as e:
                self.logger.error(f"Unexpected error in {description}: {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(self.config["parameters"]["retry_delay"])
        
        raise Exception(f"Failed after {max_retries} attempts: {description}")
    
    def setup_directories(self):
        """Create necessary directories"""
        dirs_to_create = [
            self.config["directories"]["output"],
            self.config["directories"]["temp"],
            "logs"
        ]
        
        for dir_path in dirs_to_create:
            Path(dir_path).mkdir(exist_ok=True)
            self.logger.info(f"Directory ready: {dir_path}")
    
    def build_bowtie2_index(self, transcript_filename):
        """Build bowtie2 index with error handling"""
        if "index_built" in self.progress:
            self.logger.info("Bowtie2 index already built, skipping...")
            return
            
        self.logger.info("Building Bowtie2 index...")
        command = [
            self.config["tools"]["bowtie2_build"], 
            transcript_filename, 
            "transcript"
        ]
        
        try:
            self.run_command_with_retry(command, "Building Bowtie2 index")
            self.progress["index_built"] = True
            self.save_progress()
        except Exception as e:
            self.logger.error(f"Failed to build Bowtie2 index: {str(e)}")
            raise
    
    def download_sra(self, sra_accession):
        """Download SRA file with retry logic"""
        fastq_path = Path(self.config["directories"]["temp"]) / f"{sra_accession}.fastq"
        
        if fastq_path.exists() and fastq_path.stat().st_size > 0:
            self.logger.info(f"FASTQ file already exists: {fastq_path}")
            return str(fastq_path)
        
        command = [
            self.config["tools"]["fastq_dump"], 
            sra_accession,
            "--outdir", self.config["directories"]["temp"]
        ]
        
        try:
            self.run_command_with_retry(command, f"Downloading {sra_accession}")
            
            # Verify download
            if not fastq_path.exists() or fastq_path.stat().st_size == 0:
                raise Exception(f"Downloaded file is empty or missing: {fastq_path}")
                
            return str(fastq_path)
            
        except Exception as e:
            self.logger.error(f"Failed to download {sra_accession}: {str(e)}")
            raise
    
    def align_reads(self, sra_accession, fastq_filename):
        """Align reads with bowtie2"""
        sam_filename = f"{sra_accession}_aligned_reads.sam"
        
        command = [
            self.config["tools"]["bowtie2"],
            "-x", "transcript",
            "-U", fastq_filename,
            "-S", sam_filename,
            "-p", self.config["parameters"]["bowtie2_threads"],
            "--very-sensitive-local"
        ]
        
        try:
            self.run_command_with_retry(command, f"Aligning {sra_accession}")
            return sam_filename
        except Exception as e:
            self.logger.error(f"Alignment failed for {sra_accession}: {str(e)}")
            raise
    
    def process_alignment_files(self, sra_accession, sam_filename):
        """Convert SAM to sorted BAM with indexing"""
        try:
            # Convert SAM to BAM
            bam_filename = f"{sra_accession}_aligned_reads.bam"
            with pysam.AlignmentFile(sam_filename, "r") as samfile:
                with pysam.AlignmentFile(bam_filename, "wb", template=samfile) as bamfile:
                    for read in samfile:
                        bamfile.write(read)
            
            # Sort BAM
            output_dir = Path(self.config["directories"]["output"])
            sorted_bam_filename = output_dir / f"{sra_accession}_aligned_reads.sorted.bam"
            pysam.sort("-o", str(sorted_bam_filename), bam_filename)
            
            # Index sorted BAM
            pysam.index(str(sorted_bam_filename))
            
            # Cleanup intermediate files
            Path(sam_filename).unlink(missing_ok=True)
            Path(bam_filename).unlink(missing_ok=True)
            
            # Verify final file
            if not sorted_bam_filename.exists():
                raise Exception(f"Sorted BAM file not created: {sorted_bam_filename}")
            
            # Count mapped reads for verification
            with pysam.AlignmentFile(str(sorted_bam_filename), "rb") as bamfile:
                mapped_reads = sum(1 for read in bamfile if not read.is_unmapped)
            
            self.logger.info(f"Processed {sra_accession}: {mapped_reads} mapped reads")
            return str(sorted_bam_filename), mapped_reads
            
        except Exception as e:
            self.logger.error(f"Failed to process alignment files for {sra_accession}: {str(e)}")
            raise
    
    def process_single_sra(self, sra_accession):
        """Process a single SRA file through the complete pipeline"""
        if sra_accession in self.progress["completed"]:
            self.logger.info(f"Skipping {sra_accession} (already completed)")
            return True
            
        if sra_accession in self.progress["failed"]:
            self.logger.info(f"Retrying previously failed {sra_accession}")
        
        self.logger.info(f"Processing {sra_accession}...")
        
        try:
            # Step 1: Download
            fastq_filename = self.download_sra(sra_accession)
            
            # Step 2: Align
            sam_filename = self.align_reads(sra_accession, fastq_filename)
            
            # Step 3: Convert to sorted BAM
            sorted_bam_filename, mapped_reads = self.process_alignment_files(sra_accession, sam_filename)
            
            # Step 4: Cleanup FASTQ (optional - comment out if you want to keep)
            Path(fastq_filename).unlink(missing_ok=True)
            
            # Mark as completed
            self.progress["completed"].append(sra_accession)
            if sra_accession in self.progress["failed"]:
                self.progress["failed"].remove(sra_accession)
            
            self.save_progress()
            
            self.logger.info(f"Successfully processed {sra_accession}: {mapped_reads} mapped reads")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to process {sra_accession}: {str(e)}")
            if sra_accession not in self.progress["failed"]:
                self.progress["failed"].append(sra_accession)
            self.save_progress()
            return False
    
    def run_pipeline(self, sra_accessions, transcript_filename, genes_to_process, g4_coords):
        """Run the complete pipeline"""
        start_time = datetime.now()
        self.logger.info(f"Starting pipeline with {len(sra_accessions)} SRA files")
        
        try:
            # Setup
            self.setup_directories()
            
            # Build index
            self.build_bowtie2_index(transcript_filename)
            
            # Process each SRA file
            successful = 0
            failed = 0
            
            for i, sra_accession in enumerate(sra_accessions, 1):
                self.logger.info(f"Processing {i}/{len(sra_accessions)}: {sra_accession}")
                
                if self.process_single_sra(sra_accession):
                    successful += 1
                else:
                    failed += 1
                
                # Progress update
                self.logger.info(f"Progress: {successful} successful, {failed} failed, {len(sra_accessions) - i} remaining")
            
            # Cleanup index files
            self.cleanup_index_files()
            
            # Generate analysis results if we have any successful files
            if successful > 0:
                self.logger.info("Generating analysis results...")
                self.generate_analysis_results(genes_to_process, g4_coords, sra_accessions)
            
            # Final report
            duration = datetime.now() - start_time
            self.logger.info(f"Pipeline completed in {duration}")
            self.logger.info(f"Results: {successful} successful, {failed} failed")
            
            if failed > 0:
                self.logger.warning(f"Failed files: {self.progress['failed']}")
            
        except KeyboardInterrupt:
            self.logger.info("Pipeline interrupted by user. Progress saved.")
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise
    
    def cleanup_index_files(self):
        """Clean up bowtie2 index files"""
        index_files = [f for f in os.listdir(".") if f.startswith("transcript.")]
        for index_file in index_files:
            Path(index_file).unlink(missing_ok=True)
        self.logger.info(f"Cleaned up {len(index_files)} index files")
    
    def generate_analysis_results(self, genes_to_process, g4_coords, sra_accessions):
        """Generate final analysis results"""
        try:
            all_gene_data = []
            
            for gene_name, entrez_id in genes_to_process:
                self.logger.info(f"Analyzing {gene_name} ({entrez_id})...")
                
                if entrez_id not in g4_coords:
                    self.logger.warning(f"G4 coordinates not defined for {entrez_id}. Skipping.")
                    continue
                
                try:
                    df = plotting_sorted_BAM(
                        gene_name,
                        entrez_id,
                        list_of_files=self.progress["completed"],  # Only use successful files
                        g4_coords=g4_coords
                    )
                    if df is not None and not df.empty:
                        all_gene_data.append(df)
                except Exception as e:
                    self.logger.error(f"Analysis failed for {gene_name}: {str(e)}")
            
            # Save results
            if all_gene_data:
                final_df = pd.concat(all_gene_data, ignore_index=True)
                output_file = "G4_RPBM_All_Genes.xlsx"
                final_df.to_excel(output_file, index=False)
                self.logger.info(f"Saved analysis results: {output_file}")
            else:
                self.logger.warning("No gene data was processed successfully")
                
        except Exception as e:
            self.logger.error(f"Failed to generate analysis results: {str(e)}")

def main():
    """Main function with your original parameters"""
    
    # Your original parameters
    sra_accessions = [
        'SRR970587', 'SRR970538', 'SRR970565', 'SRR970490', 'SRR970588',
        'SRR9113069', 'SRR9113063', 'SRR403885', 'SRR1248253', 'SRR1598971',
        'SRR1573934', 'SRR2064017', 'SRR2064020', 'SRR4450327'
    ]
    
    transcript_filename = "transcript11.fasta"
    
    g4_coords = {
        "NM_004565.3": (1047, 1076),
        "NM_000418.4": (2057, 2086),
        "NM_014727.3": (288, 309),
        "NM_014727.3": (3157, 3185),
        "NM_001164586.2": (7837, 7866),
        "NM_001033910.3": (1407, 1424),
        "NM_001077268.2": (614, 636),
        "NM_020832.3": (1467, 1496),
        "NM_052857.4": (522, 551),
        "NM_015355.4": (351, 380),
        "NM_024619.4": (419, 448),
        "NM_002024.6": (1810, 1838)
    }
    
    genes_to_process = [
        ("PEX14", "NM_004565.3"),
        ("IL4R", "NM_000418.4"),
        ("MLL4", "NM_014727.3"),
        ("IGFN1", "NM_001164586.2"),
        ("TRAF5", "NM_001033910.3"),
        ("ZFYVE19", "NM_001077268.2"),
        ("ZNF687", "NM_020832.3"),
        ("ZNF830", "NM_052857.4"),
        ("SUZ12", "NM_015355.4"),
        ("FN3KRP", "NM_024619.4"),
        ("FMR1", "NM_002024.6")
    ]
    
    # Run pipeline
    pipeline = RobustRNASeqPipeline()
    pipeline.run_pipeline(sra_accessions, transcript_filename, genes_to_process, g4_coords)

if __name__ == "__main__":
    main()
