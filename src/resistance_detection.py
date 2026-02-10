"""
resistance_detection.py - Drug Resistance Detection Module

This module analyzes annotated variants to classify samples as resistant
or sensitive to specific antimalarial drugs based on known resistance markers.

Author: Molecular Surveillance Team
Version: 1.0.0

Drug Resistance Markers:
- Chloroquine: pfcrt K76T, pfmdr1 N86Y
- Artemisinin: kelch13 propeller domain mutations
- SP (Sulfadoxine-Pyrimethamine): dhfr/dhps quintuple mutant
- Lumefantrine: pfmdr1 N86, 184F, D1246
- Amodiaquine: pfcrt 72-76 haplotype, pfmdr1 86Y
"""

import os
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from enum import Enum
from collections import defaultdict

from src.utils import setup_logging, ensure_directory
from src.variant_calling import KNOWN_RESISTANCE_MUTATIONS


# =============================================================================
# ENUMS AND CONSTANTS
# =============================================================================

class ResistanceStatus(Enum):
    """Resistance classification status."""
    SENSITIVE = "Sensitive"
    RESISTANT = "Resistant"
    PARTIAL = "Partial Resistance"
    UNKNOWN = "Unknown"
    
    def __str__(self):
        return self.value


class DrugClass(Enum):
    """Antimalarial drug classes."""
    AMINOQUINOLINE = "4-Aminoquinoline"
    ARTEMISININ = "Artemisinin derivatives"
    ANTIFOLATE = "Antifolate"
    ARYL_AMINO_ALCOHOL = "Aryl amino alcohol"
    
    def __str__(self):
        return self.value


# Drug classification database
DRUG_DATABASE = {
    'Chloroquine': {
        'class': DrugClass.AMINOQUINOLINE,
        'abbreviation': 'CQ',
        'markers': {
            'pfcrt': ['K76T'],  # Primary marker
            'pfmdr1': ['N86Y']  # Modulating factor
        },
        'primary_gene': 'pfcrt',
        'description': 'First-line treatment historically; widespread resistance'
    },
    'Amodiaquine': {
        'class': DrugClass.AMINOQUINOLINE,
        'abbreviation': 'AQ',
        'markers': {
            'pfcrt': ['K76T'],
            'pfmdr1': ['N86Y']
        },
        'primary_gene': 'pfcrt',
        'description': 'Used in ACT combinations (ASAQ)'
    },
    'Artemisinin': {
        'class': DrugClass.ARTEMISININ,
        'abbreviation': 'ART',
        'markers': {
            'kelch13': [
                'C580Y', 'R539T', 'Y493H', 'R561H', 'P553L',
                'I543T', 'P574L', 'N458Y', 'M476I', 'G449A',
                'Y469C', 'A675V', 'A622S'
            ]
        },
        'primary_gene': 'kelch13',
        'description': 'Core component of ACT; delayed parasite clearance'
    },
    'Lumefantrine': {
        'class': DrugClass.ARYL_AMINO_ALCOHOL,
        'abbreviation': 'LUM',
        'markers': {
            'pfmdr1': ['N86', 'Y184F', 'D1246']  # Wild-type N86 associated with reduced susceptibility
        },
        'primary_gene': 'pfmdr1',
        'description': 'Partner drug in Artemether-Lumefantrine (Coartem)'
    },
    'Mefloquine': {
        'class': DrugClass.ARYL_AMINO_ALCOHOL,
        'abbreviation': 'MQ',
        'markers': {
            'pfmdr1': ['gene_amplification', 'S1034C', 'N1042D', 'D1246Y']
        },
        'primary_gene': 'pfmdr1',
        'description': 'Used for prophylaxis and treatment'
    },
    'Pyrimethamine': {
        'class': DrugClass.ANTIFOLATE,
        'abbreviation': 'PYR',
        'markers': {
            'dhfr': ['N51I', 'C59R', 'S108N', 'I164L']
        },
        'primary_gene': 'dhfr',
        'description': 'DHFR inhibitor; component of SP'
    },
    'Sulfadoxine': {
        'class': DrugClass.ANTIFOLATE,
        'abbreviation': 'SDX',
        'markers': {
            'dhps': ['S436A', 'A437G', 'K540E', 'A581G', 'A613S']
        },
        'primary_gene': 'dhps',
        'description': 'DHPS inhibitor; component of SP'
    },
    'Sulfadoxine-Pyrimethamine': {
        'class': DrugClass.ANTIFOLATE,
        'abbreviation': 'SP',
        'markers': {
            'dhfr': ['N51I', 'C59R', 'S108N', 'I164L'],
            'dhps': ['A437G', 'K540E', 'A581G']
        },
        'primary_gene': 'dhfr',  # Both genes important
        'description': 'IPTp drug; quintuple mutant indicates resistance'
    },
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class MutationCall:
    """
    Represents a detected mutation.
    
    Attributes
    ----------
    gene : str
        Gene name
    mutation : str
        Mutation notation (e.g., 'K76T')
    position : int
        Amino acid position
    ref_aa : str
        Reference amino acid
    alt_aa : str
        Alternative amino acid
    genotype : str
        Genotype (e.g., '0/1', '1/1')
    allele_frequency : float
        Variant allele frequency
    confidence : str
        Confidence level ('High', 'Medium', 'Low')
    """
    gene: str
    mutation: str
    position: int
    ref_aa: str = ''
    alt_aa: str = ''
    genotype: str = '1/1'
    allele_frequency: float = 1.0
    confidence: str = 'High'
    
    @property
    def is_homozygous(self) -> bool:
        """Check if mutation is homozygous."""
        return self.genotype in ['1/1', '1|1']
    
    @property
    def is_heterozygous(self) -> bool:
        """Check if mutation is heterozygous."""
        return self.genotype in ['0/1', '0|1', '1/0', '1|0']


@dataclass
class DrugResistanceResult:
    """
    Resistance classification for a single drug.
    
    Attributes
    ----------
    drug : str
        Drug name
    status : ResistanceStatus
        Resistance classification
    mutations_found : list
        List of resistance mutations detected
    confidence : float
        Confidence in classification (0-1)
    notes : str
        Additional notes
    """
    drug: str
    status: ResistanceStatus = ResistanceStatus.UNKNOWN
    mutations_found: List[MutationCall] = field(default_factory=list)
    mutations_expected: List[str] = field(default_factory=list)
    confidence: float = 0.0
    notes: str = ''
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'drug': self.drug,
            'status': str(self.status),
            'mutations_found': [m.mutation for m in self.mutations_found],
            'confidence': round(self.confidence, 2),
            'notes': self.notes
        }


@dataclass
class SampleResistanceProfile:
    """
    Complete resistance profile for a sample.
    
    Attributes
    ----------
    sample_id : str
        Sample identifier
    drug_results : dict
        Results by drug
    all_mutations : list
        All resistance mutations found
    resistance_genes : set
        Genes with resistance mutations
    """
    sample_id: str
    drug_results: Dict[str, DrugResistanceResult] = field(default_factory=dict)
    all_mutations: List[MutationCall] = field(default_factory=list)
    resistance_genes: Set[str] = field(default_factory=set)
    
    @property
    def has_any_resistance(self) -> bool:
        """Check if any drug resistance detected."""
        return any(
            r.status == ResistanceStatus.RESISTANT 
            for r in self.drug_results.values()
        )
    
    @property
    def resistant_drugs(self) -> List[str]:
        """Get list of drugs with resistance."""
        return [
            drug for drug, result in self.drug_results.items()
            if result.status == ResistanceStatus.RESISTANT
        ]
    
    def get_summary(self) -> Dict:
        """Get summary of resistance profile."""
        return {
            'sample_id': self.sample_id,
            'total_mutations': len(self.all_mutations),
            'resistant_to': self.resistant_drugs,
            'genes_affected': list(self.resistance_genes),
            'drug_classifications': {
                drug: str(result.status)
                for drug, result in self.drug_results.items()
            }
        }


# =============================================================================
# RESISTANCE DETECTOR CLASS
# =============================================================================

class ResistanceDetector:
    """
    Detect and classify drug resistance from variant data.
    
    This class analyzes annotated variants to identify known resistance
    mutations and classify samples as resistant or sensitive to
    specific antimalarial drugs.
    
    Example
    -------
    >>> detector = ResistanceDetector()
    >>> profile = detector.analyze_vcf('annotated.vcf', 'sample_001')
    >>> print(f"Resistant to: {', '.join(profile.resistant_drugs)}")
    """
    
    def __init__(self):
        """Initialize resistance detector."""
        self.drug_db = DRUG_DATABASE
        self.mutation_db = KNOWN_RESISTANCE_MUTATIONS
        self.logger = setup_logging()
    
    def _parse_vcf_mutations(self, vcf_path: str) -> List[MutationCall]:
        """
        Extract resistance-relevant mutations from annotated VCF.
        
        Parameters
        ----------
        vcf_path : str
            Path to annotated VCF file
        
        Returns
        -------
        list of MutationCall
            Detected mutations
        """
        mutations = []
        
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                
                # Parse INFO field
                info = fields[7]
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                
                # Check if this is a resistance mutation
                if info_dict.get('RESISTANCE') != 'YES':
                    continue
                
                gene = info_dict.get('GENE', '')
                aa_change = info_dict.get('AA_CHANGE', '')
                
                if not gene or not aa_change:
                    continue
                
                # Parse genotype from sample column
                format_keys = fields[8].split(':')
                sample_values = fields[9].split(':')
                sample_dict = dict(zip(format_keys, sample_values))
                
                genotype = sample_dict.get('GT', '1/1')
                
                # Parse allele frequency
                af = float(info_dict.get('AF', '1.0'))
                
                # Parse mutation notation (e.g., 'K76T')
                import re
                match = re.match(r'([A-Z*])(\d+)([A-Z*])', aa_change)
                if match:
                    ref_aa = match.group(1)
                    position = int(match.group(2))
                    alt_aa = match.group(3)
                else:
                    ref_aa = ''
                    position = 0
                    alt_aa = ''
                
                mutation = MutationCall(
                    gene=gene.lower(),
                    mutation=aa_change,
                    position=position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    genotype=genotype,
                    allele_frequency=af,
                    confidence='High' if af > 0.8 else 'Medium' if af > 0.5 else 'Low'
                )
                
                mutations.append(mutation)
        
        return mutations
    
    def _classify_drug_resistance(
        self,
        drug: str,
        mutations: List[MutationCall]
    ) -> DrugResistanceResult:
        """
        Classify resistance status for a specific drug.
        
        Parameters
        ----------
        drug : str
            Drug name
        mutations : list of MutationCall
            All detected mutations
        
        Returns
        -------
        DrugResistanceResult
            Classification result
        """
        drug_info = self.drug_db.get(drug)
        if not drug_info:
            return DrugResistanceResult(
                drug=drug,
                status=ResistanceStatus.UNKNOWN,
                notes="Drug not in database"
            )
        
        result = DrugResistanceResult(drug=drug)
        
        # Get expected markers
        expected_markers = []
        for gene, markers in drug_info['markers'].items():
            for marker in markers:
                expected_markers.append((gene, marker))
        
        result.mutations_expected = [f"{g}:{m}" for g, m in expected_markers]
        
        # Find matching mutations
        found_mutations = []
        primary_gene = drug_info['primary_gene']
        primary_found = False
        
        for gene, marker in expected_markers:
            for mutation in mutations:
                if mutation.gene.lower() == gene.lower():
                    if marker in mutation.mutation:
                        found_mutations.append(mutation)
                        if gene.lower() == primary_gene.lower():
                            primary_found = True
        
        result.mutations_found = found_mutations
        
        # Determine resistance status
        if len(found_mutations) == 0:
            result.status = ResistanceStatus.SENSITIVE
            result.confidence = 0.9
            result.notes = "No resistance markers detected"
        elif primary_found:
            result.status = ResistanceStatus.RESISTANT
            result.confidence = 0.95 if len(found_mutations) > 1 else 0.85
            result.notes = f"Primary marker in {primary_gene} detected"
        else:
            result.status = ResistanceStatus.PARTIAL
            result.confidence = 0.7
            result.notes = "Secondary markers present, primary marker absent"
        
        return result
    
    def analyze_vcf(
        self,
        vcf_path: str,
        sample_id: Optional[str] = None
    ) -> SampleResistanceProfile:
        """
        Analyze VCF file for drug resistance.
        
        Parameters
        ----------
        vcf_path : str
            Path to annotated VCF file
        sample_id : str, optional
            Sample identifier
        
        Returns
        -------
        SampleResistanceProfile
            Complete resistance profile
        """
        self.logger.info(f"Analyzing resistance for: {vcf_path}")
        
        if sample_id is None:
            sample_id = os.path.splitext(os.path.basename(vcf_path))[0]
        
        # Extract mutations
        mutations = self._parse_vcf_mutations(vcf_path)
        
        # Create profile
        profile = SampleResistanceProfile(sample_id=sample_id)
        profile.all_mutations = mutations
        profile.resistance_genes = set(m.gene for m in mutations)
        
        # Classify resistance for each drug
        for drug in self.drug_db:
            result = self._classify_drug_resistance(drug, mutations)
            profile.drug_results[drug] = result
        
        self.logger.info(
            f"Analysis complete: {len(mutations)} resistance mutations, "
            f"resistant to {len(profile.resistant_drugs)} drugs"
        )
        
        return profile
    
    def generate_report(
        self,
        profile: SampleResistanceProfile,
        output_path: str
    ) -> str:
        """
        Generate a text resistance report.
        
        Parameters
        ----------
        profile : SampleResistanceProfile
            Resistance profile
        output_path : str
            Path for output report
        
        Returns
        -------
        str
            Path to generated report
        """
        lines = [
            "=" * 70,
            "ANTIMALARIAL DRUG RESISTANCE REPORT",
            "=" * 70,
            "",
            f"Sample ID: {profile.sample_id}",
            f"Total resistance mutations detected: {len(profile.all_mutations)}",
            f"Genes with mutations: {', '.join(profile.resistance_genes) or 'None'}",
            "",
            "-" * 70,
            "DRUG RESISTANCE CLASSIFICATION",
            "-" * 70,
            ""
        ]
        
        # Drug-by-drug results
        for drug, result in profile.drug_results.items():
            drug_info = self.drug_db.get(drug, {})
            
            status_symbol = {
                ResistanceStatus.SENSITIVE: "✓",
                ResistanceStatus.RESISTANT: "✗",
                ResistanceStatus.PARTIAL: "~",
                ResistanceStatus.UNKNOWN: "?"
            }.get(result.status, "?")
            
            lines.append(f"{drug} ({drug_info.get('abbreviation', '')})")
            lines.append(f"  Status: {status_symbol} {result.status}")
            lines.append(f"  Confidence: {result.confidence*100:.0f}%")
            
            if result.mutations_found:
                mutation_strs = [m.mutation for m in result.mutations_found]
                lines.append(f"  Mutations found: {', '.join(mutation_strs)}")
            else:
                lines.append("  Mutations found: None")
            
            lines.append(f"  Notes: {result.notes}")
            lines.append("")
        
        # Detailed mutation list
        lines.extend([
            "-" * 70,
            "DETAILED MUTATION LIST",
            "-" * 70,
            ""
        ])
        
        if profile.all_mutations:
            lines.append(f"{'Gene':<12} {'Mutation':<12} {'Genotype':<10} {'AF':<8} {'Confidence'}")
            lines.append("-" * 60)
            
            for mut in profile.all_mutations:
                lines.append(
                    f"{mut.gene:<12} {mut.mutation:<12} {mut.genotype:<10} "
                    f"{mut.allele_frequency:<8.2f} {mut.confidence}"
                )
        else:
            lines.append("No resistance mutations detected.")
        
        lines.extend([
            "",
            "-" * 70,
            "INTERPRETATION NOTES",
            "-" * 70,
            "",
            "• Resistant: Known resistance markers detected in primary gene",
            "• Partial: Modulating mutations present without primary marker",
            "• Sensitive: No known resistance markers detected",
            "",
            "This report is for research purposes only.",
            "Clinical decisions should incorporate additional testing.",
            "",
            "=" * 70,
            "END OF REPORT",
            "=" * 70
        ])
        
        report_text = '\n'.join(lines)
        
        ensure_directory(os.path.dirname(output_path) or '.')
        with open(output_path, 'w') as f:
            f.write(report_text)
        
        self.logger.info(f"Resistance report saved to {output_path}")
        
        return output_path


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def detect_resistance(
    vcf_path: str,
    output_dir: str,
    sample_id: Optional[str] = None
) -> Tuple[SampleResistanceProfile, str]:
    """
    Convenience function to detect resistance and generate report.
    
    Parameters
    ----------
    vcf_path : str
        Path to annotated VCF file
    output_dir : str
        Output directory for report
    sample_id : str, optional
        Sample identifier
    
    Returns
    -------
    tuple
        (SampleResistanceProfile, path to report)
    """
    ensure_directory(output_dir)
    
    if sample_id is None:
        sample_id = os.path.splitext(os.path.basename(vcf_path))[0]
        sample_id = sample_id.replace('_annotated', '')
    
    detector = ResistanceDetector()
    profile = detector.analyze_vcf(vcf_path, sample_id)
    
    report_path = os.path.join(output_dir, f"{sample_id}_resistance_report.txt")
    detector.generate_report(profile, report_path)
    
    return profile, report_path


def get_drug_info(drug_name: str) -> Optional[Dict]:
    """
    Get information about a specific drug.
    
    Parameters
    ----------
    drug_name : str
        Drug name
    
    Returns
    -------
    dict or None
        Drug information
    """
    return DRUG_DATABASE.get(drug_name)


def list_all_drugs() -> List[str]:
    """Get list of all drugs in database."""
    return list(DRUG_DATABASE.keys())


if __name__ == '__main__':
    print("Drug Resistance Detection Module")
    print("-" * 50)
    
    # Show drug database
    print("\nAntimalarial Drug Database:")
    for drug, info in DRUG_DATABASE.items():
        print(f"\n  {drug} ({info['abbreviation']})")
        print(f"    Class: {info['class']}")
        print(f"    Primary gene: {info['primary_gene']}")
        print(f"    Markers: {info['markers']}")
