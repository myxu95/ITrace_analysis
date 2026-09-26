"""
Quality Reporter

Generates comprehensive quality analysis reports and visualizations.
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union
import logging

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    logger.warning("Matplotlib/Pandas not available. Plotting features will be disabled.")


class QualityReporter:
    """
    Generates comprehensive quality analysis reports.
    
    Features:
    - Consolidated quality reports from all analysis modules
    - Statistical summaries and trends
    - Issue identification and prioritization
    - Visualization of quality metrics
    - Actionable recommendations
    """
    
    def __init__(self, output_directory: Optional[Union[str, Path]] = None):
        """
        Initialize quality reporter.
        
        Args:
            output_directory: Directory for saving reports and plots
        """
        if output_directory is None:
            self.output_dir = Path.cwd() / "quality_reports"
        else:
            self.output_dir = Path(output_directory)
            
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def generate_comprehensive_report(self, 
                                    completeness_results: List[Dict],
                                    validation_results: List[Dict],
                                    batch_tracker_data: Dict,
                                    report_name: Optional[str] = None) -> Dict:
        """
        Generate comprehensive quality report from all analysis modules.
        
        Args:
            completeness_results: Results from MDCompletenessChecker
            validation_results: Results from StructureValidator
            batch_tracker_data: Data from BatchTracker
            report_name: Optional name for the report
            
        Returns:
            Comprehensive quality report
        """
        if report_name is None:
            report_name = f"quality_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            
        logger.info(f"Generating comprehensive quality report: {report_name}")
        
        # Analyze completeness data
        completeness_summary = self._analyze_completeness_results(completeness_results)
        
        # Analyze validation data
        validation_summary = self._analyze_validation_results(validation_results)
        
        # Combine batch tracking data
        batch_summary = batch_tracker_data
        
        # Cross-analysis insights
        cross_analysis = self._perform_cross_analysis(
            completeness_results, validation_results, batch_tracker_data
        )
        
        # Generate recommendations
        recommendations = self._generate_comprehensive_recommendations(
            completeness_summary, validation_summary, batch_summary, cross_analysis
        )
        
        # Issue prioritization
        issues = self._prioritize_issues(
            completeness_summary, validation_summary, batch_summary
        )
        
        report = {
            'report_metadata': {
                'report_name': report_name,
                'generation_timestamp': datetime.now().isoformat(),
                'total_directories_analyzed': len(completeness_results),
                'total_structures_validated': len(validation_results),
                'analysis_modules': ['MDCompletenessChecker', 'StructureValidator', 'BatchTracker']
            },
            'executive_summary': self._create_executive_summary(
                completeness_summary, validation_summary, batch_summary
            ),
            'completeness_analysis': completeness_summary,
            'structure_validation': validation_summary,
            'batch_tracking': batch_summary,
            'cross_analysis': cross_analysis,
            'prioritized_issues': issues,
            'recommendations': recommendations
        }
        
        # Save report
        report_file = self.output_dir / f"{report_name}.json"
        self._save_json_report(report, report_file)
        
        # Generate visualizations if possible
        if HAS_PLOTTING:
            self._generate_visualizations(report, report_name)
        
        return report
    
    def _analyze_completeness_results(self, results: List[Dict]) -> Dict:
        """Analyze MD completeness check results."""
        if not results:
            return {'total': 0, 'analysis': 'No data available'}
        
        total = len(results)
        status_counts = {}
        score_sum = 0
        all_issues = []
        
        for result in results:
            status = result.get('status', 'unknown')
            status_counts[status] = status_counts.get(status, 0) + 1
            score_sum += result.get('completeness_score', 0)
            all_issues.extend(result.get('issues', []))
        
        issue_frequency = {}
        for issue in all_issues:
            issue_frequency[issue] = issue_frequency.get(issue, 0) + 1
        
        return {
            'total_analyzed': total,
            'status_distribution': status_counts,
            'average_completeness_score': score_sum / total if total > 0 else 0,
            'success_rate': (status_counts.get('complete', 0) / total) * 100 if total > 0 else 0,
            'common_issues': sorted(issue_frequency.items(), key=lambda x: x[1], reverse=True)[:10],
            'issue_summary': {
                'total_issues': len(all_issues),
                'unique_issue_types': len(issue_frequency),
                'most_frequent_issue': max(issue_frequency.items(), key=lambda x: x[1]) if issue_frequency else None
            }
        }
    
    def _analyze_validation_results(self, results: List[Dict]) -> Dict:
        """Analyze structure validation results."""
        if not results:
            return {'total': 0, 'analysis': 'No data available'}
        
        total = len(results)
        status_counts = {}
        score_sum = 0
        chain_counts = []
        all_issues = []
        
        for result in results:
            status = result.get('status', 'unknown')
            status_counts[status] = status_counts.get(status, 0) + 1
            score_sum += result.get('validation_score', 0)
            all_issues.extend(result.get('issues', []))
            
            # Extract chain count
            chain_analysis = result.get('chain_analysis', {})
            if 'protein_chains' in chain_analysis:
                chain_counts.append(chain_analysis['protein_chains'])
        
        issue_frequency = {}
        for issue in all_issues:
            issue_frequency[issue] = issue_frequency.get(issue, 0) + 1
        
        # Chain count analysis
        chain_stats = {}
        if chain_counts:
            chain_stats = {
                'mean_chains': sum(chain_counts) / len(chain_counts),
                'chain_distribution': {str(count): chain_counts.count(count) for count in set(chain_counts)},
                'structures_with_expected_chains': chain_counts.count(5),  # Assuming 5 is expected
                'chain_count_compliance': (chain_counts.count(5) / len(chain_counts)) * 100
            }
        
        return {
            'total_validated': total,
            'status_distribution': status_counts,
            'average_validation_score': score_sum / total if total > 0 else 0,
            'validation_rate': ((status_counts.get('valid', 0) + status_counts.get('acceptable', 0)) / total) * 100 if total > 0 else 0,
            'chain_analysis': chain_stats,
            'common_issues': sorted(issue_frequency.items(), key=lambda x: x[1], reverse=True)[:10],
            'issue_summary': {
                'total_issues': len(all_issues),
                'unique_issue_types': len(issue_frequency),
                'most_frequent_issue': max(issue_frequency.items(), key=lambda x: x[1]) if issue_frequency else None
            }
        }
    
    def _perform_cross_analysis(self, completeness_results: List[Dict], 
                              validation_results: List[Dict], 
                              batch_data: Dict) -> Dict:
        """Perform cross-module analysis to identify patterns."""
        cross_insights = {
            'correlations': {},
            'patterns': [],
            'data_quality_trends': {}
        }
        
        # Analyze score correlations if we have matching data
        if completeness_results and validation_results:
            # This would require matching by directory/file paths
            # For now, we'll analyze general trends
            
            avg_completeness = sum(r.get('completeness_score', 0) for r in completeness_results) / len(completeness_results)
            avg_validation = sum(r.get('validation_score', 0) for r in validation_results) / len(validation_results)
            
            cross_insights['correlations'] = {
                'average_completeness_score': avg_completeness,
                'average_validation_score': avg_validation,
                'quality_alignment': 'high' if abs(avg_completeness - avg_validation) < 20 else 'low'
            }
        
        # Identify patterns
        patterns = []
        
        # Pattern: High incompleteness with structure issues
        incomplete_md = len([r for r in completeness_results if r.get('status') != 'complete'])
        problematic_structures = len([r for r in validation_results if r.get('status') == 'problematic'])
        
        if incomplete_md > 0 and problematic_structures > 0:
            patterns.append({
                'pattern': 'correlated_quality_issues',
                'description': f'{incomplete_md} incomplete MD simulations and {problematic_structures} problematic structures detected',
                'severity': 'high' if (incomplete_md + problematic_structures) > len(completeness_results) * 0.3 else 'medium'
            })
        
        # Pattern: Batch processing efficiency
        if 'total_batches' in batch_data:
            avg_dirs_per_batch = batch_data.get('total_directories', 0) / max(batch_data.get('total_batches', 1), 1)
            if avg_dirs_per_batch < 50:
                patterns.append({
                    'pattern': 'small_batch_sizes',
                    'description': f'Average batch size is {avg_dirs_per_batch:.1f} directories - consider larger batches for efficiency',
                    'severity': 'low'
                })
        
        cross_insights['patterns'] = patterns
        
        return cross_insights
    
    def _generate_comprehensive_recommendations(self, completeness_summary: Dict, 
                                             validation_summary: Dict, 
                                             batch_summary: Dict, 
                                             cross_analysis: Dict) -> List[Dict]:
        """Generate actionable recommendations based on all analysis results."""
        recommendations = []
        
        # Completeness-based recommendations
        if completeness_summary.get('success_rate', 0) < 80:
            recommendations.append({
                'category': 'md_completeness',
                'priority': 'high',
                'title': 'Low MD Completion Rate',
                'description': f"Only {completeness_summary.get('success_rate', 0):.1f}% of MD simulations completed successfully",
                'actions': [
                    'Review batch job parameters and resource allocation',
                    'Check for system-wide issues or resource constraints',
                    'Implement better monitoring for long-running simulations'
                ]
            })
        
        # Validation-based recommendations
        if validation_summary.get('validation_rate', 0) < 80:
            recommendations.append({
                'category': 'structure_validation',
                'priority': 'high',
                'title': 'Structure Quality Issues',
                'description': f"Only {validation_summary.get('validation_rate', 0):.1f}% of structures passed validation",
                'actions': [
                    'Review input PDB preparation protocols',
                    'Validate starting structures before MD setup',
                    'Check for consistent system preparation parameters'
                ]
            })
        
        # Chain count issues
        chain_analysis = validation_summary.get('chain_analysis', {})
        if chain_analysis.get('chain_count_compliance', 100) < 90:
            recommendations.append({
                'category': 'structure_validation',
                'priority': 'medium',
                'title': 'Chain Count Inconsistencies',
                'description': f"Only {chain_analysis.get('chain_count_compliance', 0):.1f}% of structures have expected chain count",
                'actions': [
                    'Review PDB selection criteria',
                    'Implement pre-processing chain validation',
                    'Update expected chain count parameters if needed'
                ]
            })
        
        # Batch processing recommendations
        if batch_summary.get('cross_batch_duplicates', 0) > 0:
            recommendations.append({
                'category': 'batch_management',
                'priority': 'medium',
                'title': 'Duplicate Processing Detected',
                'description': f"{batch_summary.get('cross_batch_duplicates', 0)} PDB IDs processed multiple times across batches",
                'actions': [
                    'Implement pre-processing duplicate checking',
                    'Review batch organization strategy',
                    'Consider consolidating duplicate results'
                ]
            })
        
        # Data quality trends
        if completeness_summary.get('average_completeness_score', 0) < 70:
            recommendations.append({
                'category': 'data_quality',
                'priority': 'high',
                'title': 'Low Overall Data Quality',
                'description': f"Average completeness score is {completeness_summary.get('average_completeness_score', 0):.1f}/100",
                'actions': [
                    'Implement stricter quality gates before processing',
                    'Review and optimize MD protocols',
                    'Increase monitoring frequency during production runs'
                ]
            })
        
        return recommendations
    
    def _prioritize_issues(self, completeness_summary: Dict, 
                         validation_summary: Dict, 
                         batch_summary: Dict) -> List[Dict]:
        """Prioritize identified issues by severity and frequency."""
        issues = []
        
        # Collect and score issues
        for issue, count in completeness_summary.get('common_issues', []):
            severity = self._calculate_issue_severity(issue, count, completeness_summary.get('total_analyzed', 1))
            issues.append({
                'type': 'completeness',
                'issue': issue,
                'frequency': count,
                'severity': severity,
                'description': self._get_issue_description(issue)
            })
        
        for issue, count in validation_summary.get('common_issues', []):
            severity = self._calculate_issue_severity(issue, count, validation_summary.get('total_validated', 1))
            issues.append({
                'type': 'validation',
                'issue': issue,
                'frequency': count,
                'severity': severity,
                'description': self._get_issue_description(issue)
            })
        
        # Sort by severity and frequency
        issues.sort(key=lambda x: (x['severity'], x['frequency']), reverse=True)
        
        return issues[:15]  # Return top 15 issues
    
    def _calculate_issue_severity(self, issue: str, frequency: int, total: int) -> int:
        """Calculate issue severity score (0-100)."""
        frequency_ratio = frequency / total
        
        # Base severity by issue type
        critical_issues = ['missing_structure_file', 'missing_trajectory_file', 'no_protein_found']
        high_issues = ['simulation_not_completed', 'unreasonable_coordinates']
        medium_issues = ['trajectory_too_small', 'unusual_chain_count']
        
        base_score = 0
        if issue in critical_issues:
            base_score = 80
        elif issue in high_issues:
            base_score = 60
        elif issue in medium_issues:
            base_score = 40
        else:
            base_score = 20
        
        # Adjust by frequency
        frequency_multiplier = min(1.5, 0.5 + frequency_ratio)
        
        return int(base_score * frequency_multiplier)
    
    def _get_issue_description(self, issue: str) -> str:
        """Get human-readable description for issue codes."""
        descriptions = {
            'missing_structure_file': 'Structure file (md.gro) is missing from simulation output',
            'missing_trajectory_file': 'Trajectory file (.xtc/.trr) is missing or corrupted',
            'simulation_not_completed': 'MD simulation did not complete successfully',
            'trajectory_too_small': 'Trajectory file is smaller than expected minimum size',
            'unusual_chain_count': 'Number of protein chains differs from expected count',
            'unreasonable_coordinates': 'Atomic coordinates are outside reasonable ranges',
            'no_protein_found': 'No protein atoms detected in structure',
            'significant_missing_residues': 'Large number of residues missing from structure'
        }
        return descriptions.get(issue, f'Issue: {issue}')
    
    def _create_executive_summary(self, completeness_summary: Dict, 
                                validation_summary: Dict, 
                                batch_summary: Dict) -> Dict:
        """Create executive summary of quality analysis."""
        total_analyzed = completeness_summary.get('total_analyzed', 0) + validation_summary.get('total_validated', 0)
        
        # Overall health score (0-100)
        completeness_score = completeness_summary.get('average_completeness_score', 0)
        validation_score = validation_summary.get('average_validation_score', 0)
        overall_health = (completeness_score + validation_score) / 2
        
        # Health status
        if overall_health >= 80:
            health_status = 'excellent'
        elif overall_health >= 60:
            health_status = 'good'
        elif overall_health >= 40:
            health_status = 'fair'
        else:
            health_status = 'poor'
        
        return {
            'overall_health_score': overall_health,
            'health_status': health_status,
            'total_items_analyzed': total_analyzed,
            'successful_completion_rate': completeness_summary.get('success_rate', 0),
            'structure_validation_rate': validation_summary.get('validation_rate', 0),
            'batch_efficiency': {
                'total_batches': batch_summary.get('total_batches', 0),
                'unique_pdbs': batch_summary.get('unique_pdbs', 0),
                'duplicate_rate': (batch_summary.get('cross_batch_duplicates', 0) / max(batch_summary.get('unique_pdbs', 1), 1)) * 100
            },
            'key_concerns': self._identify_key_concerns(completeness_summary, validation_summary, batch_summary),
            'data_quality_trend': 'stable'  # Would need historical data for actual trend analysis
        }
    
    def _identify_key_concerns(self, completeness_summary: Dict, 
                             validation_summary: Dict, 
                             batch_summary: Dict) -> List[str]:
        """Identify key concerns for executive summary."""
        concerns = []
        
        if completeness_summary.get('success_rate', 0) < 70:
            concerns.append('Low MD simulation completion rate')
        
        if validation_summary.get('validation_rate', 0) < 70:
            concerns.append('High structure validation failure rate')
        
        if batch_summary.get('cross_batch_duplicates', 0) > batch_summary.get('unique_pdbs', 0) * 0.1:
            concerns.append('Significant duplicate processing detected')
        
        avg_completeness = completeness_summary.get('average_completeness_score', 0)
        avg_validation = validation_summary.get('average_validation_score', 0)
        if abs(avg_completeness - avg_validation) > 30:
            concerns.append('Inconsistent quality between completeness and validation metrics')
        
        return concerns
    
    def _save_json_report(self, report: Dict, file_path: Path):
        """Save report as JSON file."""
        try:
            with open(file_path, 'w') as f:
                json.dump(report, f, indent=2, default=str)
            logger.info(f"Quality report saved to: {file_path}")
        except Exception as e:
            logger.error(f"Failed to save report: {e}")
    
    def _generate_visualizations(self, report: Dict, report_name: str):
        """Generate visualization plots for the quality report."""
        if not HAS_PLOTTING:
            return
        
        try:
            # Create plots directory
            plots_dir = self.output_dir / f"{report_name}_plots"
            plots_dir.mkdir(exist_ok=True)
            
            # Plot 1: Status distribution pie charts
            self._plot_status_distributions(report, plots_dir)
            
            # Plot 2: Score distributions
            self._plot_score_distributions(report, plots_dir)
            
            # Plot 3: Issue frequency bar chart
            self._plot_issue_frequencies(report, plots_dir)
            
            # Plot 4: Executive summary dashboard
            self._plot_executive_dashboard(report, plots_dir)
            
            logger.info(f"Visualizations saved to: {plots_dir}")
            
        except Exception as e:
            logger.warning(f"Failed to generate visualizations: {e}")
    
    def _plot_status_distributions(self, report: Dict, plots_dir: Path):
        """Plot status distribution pie charts."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Completeness status
        completeness = report['completeness_analysis']['status_distribution']
        if completeness:
            ax1.pie(completeness.values(), labels=completeness.keys(), autopct='%1.1f%%')
            ax1.set_title('MD Completeness Status Distribution')
        
        # Validation status
        validation = report['structure_validation']['status_distribution']
        if validation:
            ax2.pie(validation.values(), labels=validation.keys(), autopct='%1.1f%%')
            ax2.set_title('Structure Validation Status Distribution')
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'status_distributions.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_score_distributions(self, report: Dict, plots_dir: Path):
        """Plot score distribution histograms."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # This would require access to individual scores
        # For now, show average scores as bar chart
        categories = ['Completeness', 'Validation']
        scores = [
            report['completeness_analysis']['average_completeness_score'],
            report['structure_validation']['average_validation_score']
        ]
        
        ax1.bar(categories, scores)
        ax1.set_ylabel('Average Score')
        ax1.set_title('Average Quality Scores')
        ax1.set_ylim(0, 100)
        
        # Health score gauge
        health_score = report['executive_summary']['overall_health_score']
        ax2.bar(['Overall Health'], [health_score], color='green' if health_score > 80 else 'orange' if health_score > 60 else 'red')
        ax2.set_ylabel('Health Score')
        ax2.set_title('Overall System Health')
        ax2.set_ylim(0, 100)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'score_distributions.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_issue_frequencies(self, report: Dict, plots_dir: Path):
        """Plot issue frequency bar chart."""
        issues = report['prioritized_issues'][:10]  # Top 10 issues
        
        if not issues:
            return
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        issue_names = [issue['issue'] for issue in issues]
        frequencies = [issue['frequency'] for issue in issues]
        severities = [issue['severity'] for issue in issues]
        
        # Color by severity
        colors = ['red' if s > 70 else 'orange' if s > 40 else 'yellow' for s in severities]
        
        bars = ax.bar(range(len(issue_names)), frequencies, color=colors)
        ax.set_xlabel('Issues')
        ax.set_ylabel('Frequency')
        ax.set_title('Top Quality Issues by Frequency')
        ax.set_xticks(range(len(issue_names)))
        ax.set_xticklabels(issue_names, rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'issue_frequencies.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_executive_dashboard(self, report: Dict, plots_dir: Path):
        """Plot executive dashboard summary."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
        
        exec_summary = report['executive_summary']
        
        # Health status
        health_score = exec_summary['overall_health_score']
        ax1.bar(['Health Score'], [health_score], color='green' if health_score > 80 else 'orange' if health_score > 60 else 'red')
        ax1.set_ylim(0, 100)
        ax1.set_title(f"Overall Health: {exec_summary['health_status'].title()}")
        
        # Success rates
        rates = ['Completion Rate', 'Validation Rate']
        values = [exec_summary['successful_completion_rate'], exec_summary['structure_validation_rate']]
        ax2.bar(rates, values)
        ax2.set_ylim(0, 100)
        ax2.set_title('Success Rates (%)')
        
        # Batch efficiency
        batch_eff = exec_summary['batch_efficiency']
        ax3.bar(['Total Batches', 'Unique PDbs'], [batch_eff['total_batches'], batch_eff['unique_pdbs']])
        ax3.set_title('Batch Processing Overview')
        
        # Key concerns
        concerns = exec_summary['key_concerns']
        if concerns:
            concern_text = '\n'.join([f"• {concern}" for concern in concerns[:5]])
            ax4.text(0.1, 0.9, 'Key Concerns:', fontweight='bold', transform=ax4.transAxes)
            ax4.text(0.1, 0.7, concern_text, transform=ax4.transAxes, verticalalignment='top')
        else:
            ax4.text(0.5, 0.5, 'No major concerns identified', ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('Key Concerns')
        ax4.axis('off')
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'executive_dashboard.png', dpi=300, bbox_inches='tight')
        plt.close()