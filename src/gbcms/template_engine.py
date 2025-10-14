"""
Thread-safe template system for output column generation.
"""

import json
import os
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Callable, Optional
from dataclasses import dataclass
from threading import Lock
import logging

logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class TemplateConfig:
    """Immutable configuration for template rendering."""
    format_name: str
    static_columns: List[str]
    count_prefix: str = "t_"
    include_strand: bool = False
    include_fragment: bool = False
    sample_columns: List[str] = None
    
    def __post_init__(self):
        if self.sample_columns is None:
            object.__setattr__(self, 'sample_columns', [])

class ColumnTemplate(ABC):
    """Abstract base class for column templates."""
    
    def __init__(self, name: str, template_type: str):
        self.name = name
        self.template_type = template_type
    
    @abstractmethod
    def generate(self, config, sample_order: List[str]) -> List[str]:
        """Generate column names based on configuration."""
        pass
    
    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}', type='{self.template_type}')"

class StaticColumnTemplate(ColumnTemplate):
    """Template for static columns that don't depend on samples or config."""
    
    def __init__(self, columns: List[str]):
        super().__init__('static_columns', 'static')
        self.columns = columns  # Immutable list
    
    def generate(self, config, sample_order: List[str]) -> List[str]:
        return self.columns.copy()

class CountColumnTemplate(ColumnTemplate):
    """Template for count-based columns with optional strand/fragment variants."""
    
    def __init__(self, prefix: str, count_types: List[str]):
        super().__init__(f'count_columns_{prefix}', 'count')
        self.prefix = prefix
        self.count_types = count_types  # Immutable list
    
    def generate(self, config, sample_order: List[str]) -> List[str]:
        columns = []
        for count_type in self.count_types:
            columns.append(f"{self.prefix}{count_type}")
            
            # Add strand-specific columns if enabled
            if config.output_strand_count:
                columns.extend([
                    f"{self.prefix}{count_type}_forward",
                    f"{self.prefix}{count_type}_reverse"
                ])
            
            # Add fragment-specific columns if enabled
            if config.output_fragment_count:
                columns.extend([
                    f"{self.prefix}{count_type}_fragment",
                    f"{self.prefix}{count_type}_fragment_forward",
                    f"{self.prefix}{count_type}_fragment_reverse"
                ])
        
        return columns

class SampleSpecificColumnTemplate(ColumnTemplate):
    """Template for columns that are specific to each sample."""
    
    def __init__(self, column_templates: List[str]):
        super().__init__('sample_columns', 'sample_specific')
        self.column_templates = column_templates  # Immutable list
    
    def generate(self, config, sample_order: List[str]) -> List[str]:
        columns = []
        for sample in sample_order:
            for template in self.column_templates:
                columns.append(f"{sample}:{template}")
                
                # Add strand-specific variants
                if config.output_strand_count:
                    columns.extend([
                        f"{sample}:{template}_FORWARD",
                        f"{sample}:{template}_REVERSE"
                    ])
                
                # Add fragment-specific variants  
                if config.output_fragment_count:
                    columns.extend([
                        f"{sample}:{template}F",  # Fragment versions
                        f"{sample}:{template}F_FORWARD",
                        f"{sample}:{template}F_REVERSE"
                    ])
        
        return columns

class TemplateEngine:
    """Thread-safe template engine for rendering column templates."""
    
    _instance = None
    _lock = Lock()
    
    def __init__(self):
        self.templates: Dict[str, List[ColumnTemplate]] = {}
        self._initialize_default_templates()
    
    @classmethod
    def get_instance(cls) -> 'TemplateEngine':
        """Get singleton instance (thread-safe)."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance
    
    def _initialize_default_templates(self):
        """Initialize default templates for all supported formats."""
        
        # VCF Template
        vcf_template = [
            StaticColumnTemplate([
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
            ]),
            SampleSpecificColumnTemplate([
                "DP", "RD", "AD", "VAF"
            ])
        ]
        
        # MAF Base Template (shared between MAF writers)
        maf_base_template = [
            StaticColumnTemplate([
                "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Variant_Classification"
            ]),
            CountColumnTemplate("t_", [
                "depth", "ref_count", "alt_count", "vaf"
            ])
        ]
        
        # Sample-Agnostic MAF Template (one row per sample per variant)
        sample_agnostic_maf_template = maf_base_template.copy()
        
        # Fillout MAF Template (all samples in one row per variant)
        fillout_maf_template = maf_base_template.copy()
        
        # Register templates
        self.register_template('vcf', vcf_template)
        self.register_template('maf', maf_base_template)
        self.register_template('sampleagnosticmaf', sample_agnostic_maf_template)
        self.register_template('fillout', fillout_maf_template)
    
    def register_template(self, format_name: str, templates: List[ColumnTemplate]):
        """Register templates for a format (thread-safe)."""
        self.templates[format_name] = templates
        logger.info(f"Registered {len(templates)} templates for format '{format_name}'")
    
    def generate_columns(self, format_name: str, config, sample_order: List[str]) -> List[str]:
        """Generate columns for a format (thread-safe, stateless)."""
        if format_name not in self.templates:
            raise ValueError(f"Unknown format: {format_name}")
        
        columns = []
        for template in self.templates[format_name]:
            template_columns = template.generate(config, sample_order)
            columns.extend(template_columns)
        
        logger.debug(f"Generated {len(columns)} columns for format '{format_name}'")
        return columns
    
    def load_template_config(self, config_file: str) -> Dict[str, Any]:
        """Load template configuration from file."""
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Template config file not found: {config_file}")
        
        with open(config_file, 'r') as f:
            return json.load(f)
    
    def get_registered_formats(self) -> List[str]:
        """Get list of registered formats (thread-safe)."""
        return list(self.templates.keys())

# Global instance for easy access

    def load_from_json_config(self, config_path: str) -> None:
        """Load template configuration from JSON file."""
        try:
            with open(config_path, r) as f:
                config_data = json.load(f)
            
            # Extract template configuration
            output_config = config_data.get(output, {})
            processing_config = config_data.get(processing, {})
            
            # Update template configurations based on JSON
            for format_name in [vcf, maf, fillout]:
                if format_name in output_config:
                    self._update_template_from_json(format_name, output_config[format_name])
                    
        except Exception as e:
            logger.error(f"Failed to load JSON config from {config_path}: {e}")
            raise

    def _update_template_from_json(self, format_name: str, format_config: dict) -> None:
        """Update template configuration from JSON format settings."""
        # This would update the template based on JSON settings
        # For now, just log the configuration
        logger.info(f"Updating {format_name} template from JSON config: {format_config}")
template_engine = TemplateEngine.get_instance()
