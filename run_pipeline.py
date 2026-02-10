#!/usr/bin/env python3
"""
Plasmodium falciparum Genomic Variant Pipeline
==============================================

Main pipeline script for molecular surveillance of antimalarial drug resistance.

This script orchestrates the complete workflow:
1. FASTQ Quality Control
2. Sequence Alignment
3. Variant Calling
4. Variant Annotation
5. Drug Resistance Detection
6. Database Storage
7. Report Generation

Usage:
    python run_pipeline.py --input sample.fastq --output results/
    python run_pipeline.py --input data/raw/ --batch --output results/

Author: Molecular Surveillance Laboratory
Version: 1.0.0
"""

import argparse
import os
import sys
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from utils import setup_logging, ensure_directory
from preprocessing import FASTQProcessor, process_fastq
from alignment import align_fastq_to_reference, ReferenceGenome
from variant_calling import call_variants_from_alignment, VariantCaller
from annotation import annotate_vcf_file, VariantAnnotator
from resistance_detection import detect_resistance, ResistanceDetector
from database import DatabaseManager, create_database


class MalariaSurveillancePipeline:
    """
    Complete pipeline for P. falciparum molecular surveillance.
    
    This pipeline processes FASTQ sequencing data through QC, alignment,
    variant calling, annotation, and resistance detection to generate
    actionable drug resistance reports for clinical and surveillance use.
    
    Attributes:
        input_path: Path to input FASTQ file or directory
        output_dir: Directory for output files
        sample_id: Sample identifier
        db_path: Path to SQLite database
        logger: Logging instance
    """
    
    def __init__(
        self,
        input_path: str,
        output_dir: str,
        sample_id: Optional[str] = None,
        db_path: Optional[str] = None,
        min_quality: int = 20,
        min_read_length: int = 50,
        verbose: bool = False
    ):
        """
        Initialize the pipeline.
        
        Args:
            input_path: Path to input FASTQ file
            output_dir: Directory for output files
            sample_id: Sample identifier (auto-generated if None)
            db_path: Path to SQLite database (default: output_dir/surveillance.db)
            min_quality: Minimum base quality for filtering
            min_read_length: Minimum read length after trimming
            verbose: Enable verbose logging
        """
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.sample_id = sample_id or self._generate_sample_id()
        self.db_path = db_path or str(self.output_dir / "surveillance.db")
        self.min_quality = min_quality
        self.min_read_length = min_read_length
        
        # Setup logging
        log_level = logging.DEBUG if verbose else logging.INFO
        self.logger = setup_logging(
            "MalariaPipeline",
            log_file=str(self.output_dir / "pipeline.log"),
            level=log_level
        )
        
        # Create output directories
        self._setup_directories()
        
        # Pipeline results storage
        self.results: Dict[str, Any] = {}
        
    def _generate_sample_id(self) -> str:
        """Generate sample ID from filename and timestamp."""
        base_name = self.input_path.stem.replace('.fastq', '').replace('.fq', '')
        timestamp = datetime.now().strftime("%Y%m%d")
        return f"{base_name}_{timestamp}"
    
    def _setup_directories(self):
        """Create output directory structure."""
        subdirs = ['qc', 'alignment', 'variants', 'annotation', 'reports']
        for subdir in subdirs:
            ensure_directory(self.output_dir / subdir)
        self.logger.info(f"Output directories created: {self.output_dir}")
    
    def run(self) -> Dict[str, Any]:
        """
        Execute the complete pipeline.
        
        Returns:
            Dictionary containing all pipeline results and paths
        """
        self.logger.info("=" * 60)
        self.logger.info("PLASMODIUM FALCIPARUM MOLECULAR SURVEILLANCE PIPELINE")
        self.logger.info("=" * 60)
        self.logger.info(f"Sample ID: {self.sample_id}")
        self.logger.info(f"Input: {self.input_path}")
        self.logger.info(f"Output: {self.output_dir}")
        self.logger.info("=" * 60)
        
        start_time = datetime.now()
        
        try:
            # Step 1: Quality Control
            self._run_qc()
            
            # Step 2: Alignment
            self._run_alignment()
            
            # Step 3: Variant Calling
            self._run_variant_calling()
            
            # Step 4: Annotation
            self._run_annotation()
            
            # Step 5: Resistance Detection
            self._run_resistance_detection()
            
            # Step 6: Database Storage
            self._store_to_database()
            
            # Step 7: Generate Reports
            self._generate_reports()
            
            # Calculate runtime
            end_time = datetime.now()
            runtime = (end_time - start_time).total_seconds()
            
            self.results['runtime_seconds'] = runtime
            self.results['status'] = 'SUCCESS'
            
            self.logger.info("=" * 60)
            self.logger.info(f"PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info(f"Total runtime: {runtime:.2f} seconds")
            self.logger.info("=" * 60)
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            self.results['status'] = 'FAILED'
            self.results['error'] = str(e)
            raise
        
        return self.results
    
    def _run_qc(self):
        """Step 1: FASTQ Quality Control."""
        self.logger.info("\n[STEP 1/7] FASTQ QUALITY CONTROL")
        self.logger.info("-" * 40)
        
        # Process FASTQ
        processor = FASTQProcessor(str(self.input_path))
        
        # Calculate QC statistics
        qc_stats = processor.get_quality_stats()
        
        self.logger.info(f"Total reads: {qc_stats.total_reads:,}")
        self.logger.info(f"Total bases: {qc_stats.total_bases:,}")
        self.logger.info(f"Mean quality: {qc_stats.mean_quality:.2f}")
        self.logger.info(f"GC content: {qc_stats.gc_content:.1f}%")
        self.logger.info(f"Q20 bases: {qc_stats.q20_percentage:.1f}%")
        self.logger.info(f"Q30 bases: {qc_stats.q30_percentage:.1f}%")
        
        # Filter reads
        filtered_path = self.output_dir / "qc" / f"{self.sample_id}_filtered.fastq"
        passed, failed = processor.filter_reads(
            output_path=str(filtered_path),
            min_quality=self.min_quality,
            min_length=self.min_read_length
        )
        
        self.logger.info(f"Reads passed filter: {passed:,} ({100*passed/(passed+failed):.1f}%)")
        self.logger.info(f"Reads failed filter: {failed:,}")
        
        # Generate QC report
        qc_report_path = self.output_dir / "qc" / f"{self.sample_id}_qc_report.txt"
        processor.generate_report(str(qc_report_path))
        
        # Store results
        self.results['qc'] = {
            'stats': qc_stats,
            'reads_passed': passed,
            'reads_failed': failed,
            'filtered_fastq': str(filtered_path),
            'report': str(qc_report_path)
        }
        
        self.logger.info(f"QC report saved: {qc_report_path}")
    
    def _run_alignment(self):
        """Step 2: Sequence Alignment."""
        self.logger.info("\n[STEP 2/7] SEQUENCE ALIGNMENT")
        self.logger.info("-" * 40)
        
        filtered_fastq = self.results['qc']['filtered_fastq']
        sam_path = self.output_dir / "alignment" / f"{self.sample_id}.sam"
        bam_path = self.output_dir / "alignment" / f"{self.sample_id}.bam"
        
        # Run alignment
        alignment_stats = align_fastq_to_reference(
            fastq_path=filtered_fastq,
            output_sam=str(sam_path),
            output_bam=str(bam_path)
        )
        
        self.logger.info(f"Total reads aligned: {alignment_stats.total_reads:,}")
        self.logger.info(f"Mapped reads: {alignment_stats.mapped_reads:,} ({alignment_stats.mapping_rate:.1f}%)")
        self.logger.info(f"Properly paired: {alignment_stats.properly_paired:,}")
        self.logger.info(f"Mean mapping quality: {alignment_stats.mean_mapping_quality:.1f}")
        
        # Store results
        self.results['alignment'] = {
            'stats': alignment_stats,
            'sam_path': str(sam_path),
            'bam_path': str(bam_path)
        }
        
        self.logger.info(f"Alignment saved: {bam_path}")
    
    def _run_variant_calling(self):
        """Step 3: Variant Calling."""
        self.logger.info("\n[STEP 3/7] VARIANT CALLING")
        self.logger.info("-" * 40)
        
        sam_path = self.results['alignment']['sam_path']
        vcf_path = self.output_dir / "variants" / f"{self.sample_id}.vcf"
        filtered_vcf_path = self.output_dir / "variants" / f"{self.sample_id}_filtered.vcf"
        
        # Call variants
        caller = VariantCaller()
        variant_stats = caller.call_variants(
            sam_path=sam_path,
            output_vcf=str(vcf_path),
            sample_id=self.sample_id
        )
        
        self.logger.info(f"Total variants called: {variant_stats.total_variants:,}")
        self.logger.info(f"SNPs: {variant_stats.snps:,}")
        self.logger.info(f"Insertions: {variant_stats.insertions:,}")
        self.logger.info(f"Deletions: {variant_stats.deletions:,}")
        
        # Filter variants
        caller.filter_vcf(
            input_vcf=str(vcf_path),
            output_vcf=str(filtered_vcf_path),
            min_quality=20,
            min_depth=10
        )
        
        # Store results
        self.results['variant_calling'] = {
            'stats': variant_stats,
            'vcf_path': str(vcf_path),
            'filtered_vcf_path': str(filtered_vcf_path)
        }
        
        self.logger.info(f"VCF saved: {filtered_vcf_path}")
    
    def _run_annotation(self):
        """Step 4: Variant Annotation."""
        self.logger.info("\n[STEP 4/7] VARIANT ANNOTATION")
        self.logger.info("-" * 40)
        
        vcf_path = self.results['variant_calling']['filtered_vcf_path']
        annotated_vcf_path = self.output_dir / "annotation" / f"{self.sample_id}_annotated.vcf"
        
        # Annotate variants
        annotator = VariantAnnotator()
        annotation_stats = annotator.annotate_vcf(
            input_vcf=vcf_path,
            output_vcf=str(annotated_vcf_path)
        )
        
        self.logger.info(f"Variants annotated: {annotation_stats.total_annotated:,}")
        self.logger.info(f"Missense mutations: {annotation_stats.missense:,}")
        self.logger.info(f"Synonymous mutations: {annotation_stats.synonymous:,}")
        self.logger.info(f"Resistance mutations: {annotation_stats.resistance_mutations:,}")
        
        # Get gene summary
        gene_summary = annotator.get_gene_summary()
        for gene, counts in gene_summary.items():
            self.logger.info(f"  {gene}: {counts['total']} variants ({counts['resistance']} resistance)")
        
        # Store results
        self.results['annotation'] = {
            'stats': annotation_stats,
            'annotated_vcf_path': str(annotated_vcf_path),
            'gene_summary': gene_summary
        }
        
        self.logger.info(f"Annotated VCF saved: {annotated_vcf_path}")
    
    def _run_resistance_detection(self):
        """Step 5: Drug Resistance Detection."""
        self.logger.info("\n[STEP 5/7] DRUG RESISTANCE DETECTION")
        self.logger.info("-" * 40)
        
        annotated_vcf = self.results['annotation']['annotated_vcf_path']
        
        # Detect resistance
        detector = ResistanceDetector()
        resistance_profile = detector.analyze_vcf(
            vcf_path=annotated_vcf,
            sample_id=self.sample_id
        )
        
        self.logger.info(f"Mutations detected: {len(resistance_profile.mutations_detected)}")
        
        # Log resistance results for each drug
        for result in resistance_profile.drug_results:
            status_str = result.status.value.upper()
            if result.status.value == 'resistant':
                self.logger.warning(f"  {result.drug}: {status_str}")
            else:
                self.logger.info(f"  {result.drug}: {status_str}")
        
        # Generate resistance report
        report_path = self.output_dir / "reports" / f"{self.sample_id}_resistance_report.txt"
        report = detector.generate_report(resistance_profile)
        with open(report_path, 'w') as f:
            f.write(report)
        
        # Store results
        self.results['resistance'] = {
            'profile': resistance_profile,
            'report_path': str(report_path),
            'has_resistance': resistance_profile.has_any_resistance,
            'resistant_drugs': resistance_profile.resistant_drugs
        }
        
        self.logger.info(f"Resistance report saved: {report_path}")
        
        if resistance_profile.has_any_resistance:
            self.logger.warning(f"⚠️  RESISTANCE DETECTED: {', '.join(resistance_profile.resistant_drugs)}")
    
    def _store_to_database(self):
        """Step 6: Store Results in Database."""
        self.logger.info("\n[STEP 6/7] DATABASE STORAGE")
        self.logger.info("-" * 40)
        
        # Initialize database
        db = create_database(self.db_path)
        
        # Insert sample record
        qc_stats = self.results['qc']['stats']
        align_stats = self.results['alignment']['stats']
        
        sample_data = {
            'sample_id': self.sample_id,
            'collection_date': datetime.now().strftime("%Y-%m-%d"),
            'location': 'Laboratory',  # Can be parameterized
            'country': 'Unknown',  # Can be parameterized
            'sequencing_platform': 'Illumina',
            'total_reads': qc_stats.total_reads,
            'mean_quality': qc_stats.mean_quality,
            'coverage': align_stats.mapping_rate
        }
        
        db.insert_sample(**sample_data)
        self.logger.info(f"Sample record inserted: {self.sample_id}")
        
        # Import variants from VCF
        annotated_vcf = self.results['annotation']['annotated_vcf_path']
        variant_count = db.import_vcf(annotated_vcf, self.sample_id)
        self.logger.info(f"Variants imported: {variant_count}")
        
        # Store resistance profile
        resistance_profile = self.results['resistance']['profile']
        db.store_resistance_profile(resistance_profile)
        self.logger.info(f"Resistance calls stored: {len(resistance_profile.drug_results)}")
        
        # Store results
        self.results['database'] = {
            'db_path': self.db_path,
            'variants_stored': variant_count
        }
        
        self.logger.info(f"Database updated: {self.db_path}")
    
    def _generate_reports(self):
        """Step 7: Generate Summary Reports."""
        self.logger.info("\n[STEP 7/7] REPORT GENERATION")
        self.logger.info("-" * 40)
        
        # Generate summary report
        summary_path = self.output_dir / "reports" / f"{self.sample_id}_summary.txt"
        
        with open(summary_path, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("PLASMODIUM FALCIPARUM MOLECULAR SURVEILLANCE REPORT\n")
            f.write("=" * 70 + "\n\n")
            
            # Sample Information
            f.write("SAMPLE INFORMATION\n")
            f.write("-" * 40 + "\n")
            f.write(f"Sample ID:        {self.sample_id}\n")
            f.write(f"Analysis Date:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input File:       {self.input_path.name}\n\n")
            
            # QC Summary
            qc = self.results['qc']
            f.write("SEQUENCING QUALITY\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Reads:      {qc['stats'].total_reads:,}\n")
            f.write(f"Reads Passed QC:  {qc['reads_passed']:,}\n")
            f.write(f"Mean Quality:     {qc['stats'].mean_quality:.1f}\n")
            f.write(f"GC Content:       {qc['stats'].gc_content:.1f}%\n\n")
            
            # Alignment Summary
            align = self.results['alignment']
            f.write("ALIGNMENT STATISTICS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Mapping Rate:     {align['stats'].mapping_rate:.1f}%\n")
            f.write(f"Mapped Reads:     {align['stats'].mapped_reads:,}\n\n")
            
            # Variant Summary
            var = self.results['variant_calling']
            f.write("VARIANT CALLING\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Variants:   {var['stats'].total_variants}\n")
            f.write(f"SNPs:             {var['stats'].snps}\n")
            f.write(f"Indels:           {var['stats'].insertions + var['stats'].deletions}\n\n")
            
            # Annotation Summary
            annot = self.results['annotation']
            f.write("VARIANT ANNOTATION\n")
            f.write("-" * 40 + "\n")
            f.write(f"Annotated:        {annot['stats'].total_annotated}\n")
            f.write(f"Missense:         {annot['stats'].missense}\n")
            f.write(f"Resistance:       {annot['stats'].resistance_mutations}\n\n")
            
            # Resistance Summary
            res = self.results['resistance']
            f.write("DRUG RESISTANCE PROFILE\n")
            f.write("-" * 40 + "\n")
            
            if res['has_resistance']:
                f.write("⚠️  RESISTANCE DETECTED\n\n")
                f.write(f"Resistant to: {', '.join(res['resistant_drugs'])}\n\n")
            else:
                f.write("✓ No resistance markers detected\n\n")
            
            f.write("Drug-by-Drug Results:\n")
            for result in res['profile'].drug_results:
                status = result.status.value.upper()
                marker = "⚠️ " if result.status.value == 'resistant' else "  "
                f.write(f"{marker}{result.drug:20s} {status}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("END OF REPORT\n")
            f.write("=" * 70 + "\n")
        
        self.results['reports'] = {
            'summary_path': str(summary_path)
        }
        
        self.logger.info(f"Summary report saved: {summary_path}")
        
        # Print final summary to console
        print("\n" + "=" * 60)
        print("ANALYSIS COMPLETE")
        print("=" * 60)
        print(f"Sample: {self.sample_id}")
        print(f"Variants: {self.results['variant_calling']['stats'].total_variants}")
        print(f"Resistance mutations: {self.results['annotation']['stats'].resistance_mutations}")
        
        if self.results['resistance']['has_resistance']:
            print(f"\n⚠️  RESISTANCE DETECTED: {', '.join(self.results['resistance']['resistant_drugs'])}")
        else:
            print("\n✓ No drug resistance detected")
        
        print(f"\nReports: {self.output_dir / 'reports'}")
        print("=" * 60)


def generate_synthetic_data(output_path: str, num_reads: int = 1000) -> str:
    """
    Generate synthetic FASTQ data for testing.
    
    Creates realistic P. falciparum-like sequencing reads for pipeline testing.
    
    Args:
        output_path: Path for output FASTQ file
        num_reads: Number of reads to generate
        
    Returns:
        Path to generated FASTQ file
    """
    import random
    
    # P. falciparum has ~19% GC content
    nucleotides = 'ATATATAT' + 'GC'  # Biased toward AT
    
    # Drug resistance gene regions to simulate reads from
    target_regions = [
        ('pfcrt', 'Pf3D7_07_v3', 403000, 406000),
        ('pfmdr1', 'Pf3D7_05_v3', 957000, 962000),
        ('kelch13', 'Pf3D7_13_v3', 1724000, 1726000),
        ('dhfr', 'Pf3D7_04_v3', 748000, 750000),
        ('dhps', 'Pf3D7_08_v3', 549000, 551000),
    ]
    
    with open(output_path, 'w') as f:
        for i in range(num_reads):
            # Generate read ID
            read_id = f"@SIM:{i+1:06d}:FLOWCELL:1:1:1:{i+1}"
            
            # Randomly select target region
            gene, chrom, start, end = random.choice(target_regions)
            pos = random.randint(start, end)
            
            # Generate sequence (100-150 bp)
            read_length = random.randint(100, 150)
            sequence = ''.join(random.choice(nucleotides) for _ in range(read_length))
            
            # Add some Ns randomly (0-2%)
            seq_list = list(sequence)
            for j in range(len(seq_list)):
                if random.random() < 0.01:
                    seq_list[j] = 'N'
            sequence = ''.join(seq_list)
            
            # Generate quality scores (Phred+33)
            # Simulate decreasing quality at read ends
            qualities = []
            for j in range(read_length):
                if j < 10:  # First 10 bp slightly lower
                    q = random.randint(25, 38)
                elif j > read_length - 15:  # Last 15 bp lower
                    q = random.randint(20, 35)
                else:
                    q = random.randint(30, 40)
                qualities.append(chr(q + 33))
            quality_string = ''.join(qualities)
            
            # Add region info to header
            header = f"{read_id} {chrom}:{pos} length={read_length}"
            
            # Write FASTQ record
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")
            f.write("+\n")
            f.write(f"{quality_string}\n")
    
    return output_path


def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Plasmodium falciparum Molecular Surveillance Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on single sample
  python run_pipeline.py -i sample.fastq -o results/
  
  # Generate synthetic test data and run
  python run_pipeline.py --generate-test-data -o results/
  
  # Run with custom parameters
  python run_pipeline.py -i sample.fastq -o results/ -q 25 -l 75 -v

For more information, visit: https://github.com/mykiie85/Genomic-Variant-Analysis-Pipeline-for-Drug-Resistant-Malaria-Surveillance
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        help='Input FASTQ file'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='./output',
        help='Output directory (default: ./output)'
    )
    
    parser.add_argument(
        '-s', '--sample-id',
        type=str,
        help='Sample identifier (auto-generated if not provided)'
    )
    
    parser.add_argument(
        '-d', '--database',
        type=str,
        help='Path to SQLite database'
    )
    
    parser.add_argument(
        '-q', '--min-quality',
        type=int,
        default=20,
        help='Minimum base quality for filtering (default: 20)'
    )
    
    parser.add_argument(
        '-l', '--min-length',
        type=int,
        default=50,
        help='Minimum read length (default: 50)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--generate-test-data',
        action='store_true',
        help='Generate synthetic test data and run pipeline'
    )
    
    parser.add_argument(
        '--num-reads',
        type=int,
        default=1000,
        help='Number of reads for synthetic data (default: 1000)'
    )
    
    args = parser.parse_args()
    
    # Handle test data generation
    if args.generate_test_data:
        ensure_directory(args.output)
        test_fastq = os.path.join(args.output, 'synthetic_sample.fastq')
        print(f"Generating synthetic test data: {test_fastq}")
        generate_synthetic_data(test_fastq, args.num_reads)
        args.input = test_fastq
        args.sample_id = args.sample_id or 'SYNTHETIC_TEST_001'
    
    # Validate input
    if not args.input:
        parser.error("Input file required. Use -i/--input or --generate-test-data")
    
    if not os.path.exists(args.input):
        parser.error(f"Input file not found: {args.input}")
    
    # Run pipeline
    pipeline = MalariaSurveillancePipeline(
        input_path=args.input,
        output_dir=args.output,
        sample_id=args.sample_id,
        db_path=args.database,
        min_quality=args.min_quality,
        min_read_length=args.min_length,
        verbose=args.verbose
    )
    
    results = pipeline.run()
    
    return 0 if results['status'] == 'SUCCESS' else 1


if __name__ == '__main__':
    sys.exit(main())
