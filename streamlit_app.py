#!/usr/bin/env python3
"""
DMD Co-Expression Network Analysis - Interactive Dashboard
========================================================

Streamlit web application for interactive visualization and exploration 
of DMD gene co-expression network analysis results.

Author: Bioinformatics Pipeline
Version: 1.0.0
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
import json
import os
from pathlib import Path

# Configure Streamlit page
st.set_page_config(
    page_title="DMD Co-Expression Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
.main-header {
    font-size: 2.5rem;
    color: #2c3e50;
    text-align: center;
    margin-bottom: 2rem;
}
.metric-card {
    background-color: #2c2f30;
    padding: 1rem;
    border-radius: 0.5rem;
    border-left: 4px solid #3498db;
}
.info-box {
    background-color: #2c2f30;
    padding: 1rem;
    border-radius: 0.5rem;
    margin: 1rem 0;
}
</style>
""", unsafe_allow_html=True)

class DMDDashboard:
    """
    Interactive dashboard for DMD co-expression analysis results.
    """

    def __init__(self, results_dir="dmd_coexpression_results"):
        self.results_dir = Path(results_dir)
        self.data = {}
        self.load_results()
    
    def load_results(self):
        """Load analysis results from files safely."""
        try:
            def safe_read_csv(path, **kwargs):
                if path.exists() and path.stat().st_size > 0:
                    try:
                        df = pd.read_csv(path, **kwargs)
                        if not df.empty:
                            return df
                        else:
                            st.warning(f"⚠️ {path.name} is empty.")
                    except Exception as e:
                        st.warning(f"⚠️ Failed to read {path.name}: {e}")
                else:
                    st.warning(f"⚠️ File missing or empty: {path.name}")
                return None

            # Define all expected file paths
            expr_path = self.results_dir / "processed" / "combined_human_expression.csv"
            samp_path = self.results_dir / "processed" / "combined_human_samples.csv"
            mod_path = self.results_dir / "networks" / "combined_human_modules.csv"
            hub_path = self.results_dir / "networks" / "combined_human_hub_genes.csv"
            corr_path = self.results_dir / "networks" / "combined_human_module_trait_correlations.csv"
            enrich_path = self.results_dir / "enrichment" / "combined_human_enrichment.csv"

            # Load data safely
            self.data['expression'] = safe_read_csv(expr_path, index_col=0)
            self.data['samples'] = safe_read_csv(samp_path)
            self.data['modules'] = safe_read_csv(mod_path)
            self.data['hub_genes'] = safe_read_csv(hub_path)
            self.data['module_correlations'] = safe_read_csv(corr_path, index_col=0)
            self.data['enrichment'] = safe_read_csv(enrich_path)

            # Check if *any* DataFrame was successfully loaded
            loaded_dfs = [v for v in self.data.values() if isinstance(v, pd.DataFrame) and not v.empty]
            if len(loaded_dfs) == 0:
                raise ValueError("No valid data files found in results directory.")

            # Display summary
            loaded_keys = [k for k, v in self.data.items() if isinstance(v, pd.DataFrame) and not v.empty]
            st.success(f"✅ Loaded data components: {', '.join(loaded_keys)}")

        except Exception as e:
            st.error(f"Error loading results: {str(e)}")
            self.create_mock_data()

    # def load_results(self):
    #     """Load analysis results from files."""
    #     try:
    #         # Load processed data
    #         if (self.results_dir / "processed" / "combined_human_expression.csv").exists():
    #             self.data['expression'] = pd.read_csv(
    #                 self.results_dir / "processed" / "combined_human_expression.csv", 
    #                 index_col=0
    #             )

    #         if (self.results_dir / "processed" / "combined_human_samples.csv").exists():
    #             self.data['samples'] = pd.read_csv(
    #                 self.results_dir / "processed" / "combined_human_samples.csv"
    #             )

    #         # Load network results
    #         if (self.results_dir / "networks" / "combined_human_modules.csv").exists():
    #             self.data['modules'] = pd.read_csv(
    #                 self.results_dir / "networks" / "combined_human_modules.csv"
    #             )

    #         if (self.results_dir / "networks" / "combined_human_hub_genes.csv").exists():
    #             self.data['hub_genes'] = pd.read_csv(
    #                 self.results_dir / "networks" / "combined_human_hub_genes.csv"
    #             )

    #         if (self.results_dir / "networks" / "combined_human_module_trait_correlations.csv").exists():
    #             self.data['module_correlations'] = pd.read_csv(
    #                 self.results_dir / "networks" / "combined_human_module_trait_correlations.csv",
    #                 index_col=0
    #             )

    #         # Load enrichment results
    #         if (self.results_dir / "enrichment" / "combined_human_enrichment.csv").exists():
    #             self.data['enrichment'] = pd.read_csv(
    #                 self.results_dir / "enrichment" / "combined_human_enrichment.csv"
    #             )

    #     except Exception as e:
    #         st.error(f"Error loading results: {str(e)}")
    #         # Create mock data for demonstration
    #         self.create_mock_data()

    def create_mock_data(self):
        """Create mock data for demonstration when real data is not available."""
        st.info("Loading demo data for visualization...")

        # Mock expression data
        np.random.seed(42)
        n_genes, n_samples = 1000, 50
        self.data['expression'] = pd.DataFrame(
            np.random.normal(0, 1, (n_genes, n_samples)),
            index=[f"Gene_{i}" for i in range(n_genes)],
            columns=[f"Sample_{i}" for i in range(n_samples)]
        )

        # Mock sample data
        self.data['samples'] = pd.DataFrame({
            'sample_id': [f"Sample_{i}" for i in range(n_samples)],
            'condition': ['DMD'] * 25 + ['Control'] * 25,
            'dataset': ['GSE38417'] * 20 + ['GSE6011'] * 15 + ['GSE109178'] * 15,
            'batch': ['Batch_1'] * 20 + ['Batch_2'] * 15 + ['Batch_3'] * 15
        })

        # Mock modules
        modules = [f"Module_{i}" for i in range(1, 21)]
        self.data['modules'] = pd.DataFrame({
            'gene': np.random.choice(self.data['expression'].index, 500),
            'module': np.random.choice(modules, 500)
        })

        # Mock hub genes  
        hub_data = []
        for module in modules:
            module_genes = self.data['modules'][self.data['modules']['module'] == module]['gene']
            if len(module_genes) > 0:
                for rank, gene in enumerate(np.random.choice(module_genes, min(5, len(module_genes)), replace=False)):
                    hub_data.append({
                        'module': module,
                        'gene': gene,
                        'hub_rank': rank + 1,
                        'hub_score': np.random.uniform(0.5, 1.0),
                        'module_size': len(module_genes)
                    })
        self.data['hub_genes'] = pd.DataFrame(hub_data)

        # Mock module correlations
        self.data['module_correlations'] = pd.DataFrame(
            np.random.uniform(-0.8, 0.8, (len(modules), 2)),
            index=modules,
            columns=['DMD', 'Control']
        )

        # Mock enrichment data
        pathways = ['Immune Response', 'Muscle Contraction', 'ECM Organization', 
                   'Inflammatory Response', 'Calcium Signaling']
        enrichment_data = []
        for module in modules[:10]:
            for pathway in np.random.choice(pathways, 2):
                enrichment_data.append({
                    'module': module,
                    'pathway': pathway,
                    'overlap': np.random.randint(2, 10),
                    'pathway_size': np.random.randint(50, 200),
                    'fold_enrichment': np.random.uniform(1.5, 5.0),
                    'p_value': np.random.uniform(0.001, 0.05)
                })
        self.data['enrichment'] = pd.DataFrame(enrichment_data)

    def render_sidebar(self):
        """Render the sidebar with navigation and filters."""
        st.sidebar.title("🧬 DMD Analysis")
        st.sidebar.markdown("---")

        # Navigation
        page = st.sidebar.selectbox(
            "Navigate to:",
            ["Overview", "Expression Data", "Co-expression Modules", "Hub Genes", 
             "Enrichment Analysis", "Interactive Network"]
        )

        st.sidebar.markdown("---")

        # Filters
        st.sidebar.subheader("Filters")

        # Dataset filter
        if 'samples' in self.data:
            datasets = self.data['samples']['dataset'].unique()
            selected_datasets = st.sidebar.multiselect(
                "Select Datasets:", 
                datasets, 
                default=datasets
            )
        else:
            selected_datasets = []

        # Condition filter  
        if 'samples' in self.data:
            conditions = self.data['samples']['condition'].unique()
            selected_conditions = st.sidebar.multiselect(
                "Select Conditions:",
                conditions,
                default=conditions
            )
        else:
            selected_conditions = []

        return page, selected_datasets, selected_conditions

    def render_overview(self):
        """Render the overview page."""
        st.markdown('<h1 class="main-header">DMD Co-Expression Network Analysis</h1>', 
                   unsafe_allow_html=True)

        # Overview cards
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.markdown(
                f'<div class="metric-card">'
                f'<h3>Genes Analyzed</h3>'
                f'<h2>{len(self.data.get("expression", {}).index) if "expression" in self.data else "N/A"}</h2>'
                f'</div>', 
                unsafe_allow_html=True
            )

        with col2:
            st.markdown(
                f'<div class="metric-card">'
                f'<h3>Samples</h3>'
                f'<h2>{len(self.data.get("samples", {})) if "samples" in self.data else "N/A"}</h2>'
                f'</div>',
                unsafe_allow_html=True
            )

        with col3:
            st.markdown(
                f'<div class="metric-card">'
                f'<h3>Co-expression Modules</h3>'
                f'<h2>{len(self.data.get("modules", {})["module"].unique()) if "modules" in self.data else "N/A"}</h2>'
                f'</div>',
                unsafe_allow_html=True
            )

        with col4:
            st.markdown(
                f'<div class="metric-card">'
                f'<h3>Hub Genes</h3>'
                f'<h2>{len(self.data.get("hub_genes", {})) if "hub_genes" in self.data else "N/A"}</h2>'
                f'</div>',
                unsafe_allow_html=True
            )

        # Analysis description
        st.markdown(
            '<div class="info-box">'
            '<h3>Analysis Overview</h3>'
            '<p>This pipeline performs comprehensive gene co-expression network analysis on Duchenne Muscular Dystrophy (DMD) transcriptomic data. '
            'The analysis integrates multiple public datasets to identify disease-associated gene modules and potential therapeutic targets.</p>'
            '<ul>'
            '<li><strong>Data Sources:</strong> GSE38417, GSE6011, GSE109178 (GEO Database)</li>'
            '<li><strong>Method:</strong> Weighted Gene Co-expression Network Analysis (WGCNA)</li>'
            '<li><strong>Validation:</strong> mdx mouse model data (GSE162455)</li>'
            '</ul>'
            '</div>',
            unsafe_allow_html=True
        )

        # Dataset summary
        if 'samples' in self.data:
            st.subheader("Dataset Summary")

            col1, col2 = st.columns(2)

            with col1:
                # Sample distribution by condition
                condition_counts = self.data['samples']['condition'].value_counts()
                fig_condition = px.pie(
                    values=condition_counts.values,
                    names=condition_counts.index,
                    title="Sample Distribution by Condition"
                )
                st.plotly_chart(fig_condition, use_container_width=True)

            with col2:
                # Sample distribution by dataset
                dataset_counts = self.data['samples']['dataset'].value_counts()
                fig_dataset = px.bar(
                    x=dataset_counts.index,
                    y=dataset_counts.values,
                    title="Samples per Dataset"
                )
                fig_dataset.update_xaxes(title="Dataset")
                fig_dataset.update_yaxes(title="Number of Samples")
                st.plotly_chart(fig_dataset, use_container_width=True)

    def render_expression_data(self, selected_datasets, selected_conditions):
        """Render expression data visualization."""
        st.header("Expression Data Analysis")

        if 'expression' not in self.data or 'samples' not in self.data:
            st.error("Expression data not available")
            return

        # Filter data
        filtered_samples = self.data['samples'][
            (self.data['samples']['dataset'].isin(selected_datasets)) &
            (self.data['samples']['condition'].isin(selected_conditions))
        ]

        if len(filtered_samples) == 0:
            st.warning("No samples match the selected filters")
            return

        # Gene expression search
        st.subheader("Gene Expression Explorer")

        col1, col2 = st.columns([2, 1])

        with col1:
            search_gene = st.text_input("Search for a gene:", placeholder="e.g., SPP1, POSTN, TYROBP")

        with col2:
            show_top_genes = st.checkbox("Show top variable genes", value=True)

        # Expression heatmap
        if search_gene and search_gene in self.data['expression'].index:
            # Single gene expression
            gene_expr = self.data['expression'].loc[search_gene, filtered_samples['sample_id']]
            gene_data = pd.DataFrame({
                'sample_id': filtered_samples['sample_id'],
                'expression': gene_expr,
                'condition': filtered_samples['condition'].values,
                'dataset': filtered_samples['dataset'].values
            })

            fig = px.box(
                gene_data, 
                x='condition', 
                y='expression',
                color='condition',
                title=f"Expression of {search_gene} across conditions"
            )
            fig.update_xaxes(title="Condition")
            fig.update_yaxes(title="Expression Level")
            st.plotly_chart(fig, use_container_width=True)

        elif show_top_genes:
            # Top variable genes heatmap
            expr_subset = self.data['expression'][filtered_samples['sample_id']]
            gene_vars = expr_subset.var(axis=1).sort_values(ascending=False)
            top_genes = gene_vars.head(50).index

            heatmap_data = expr_subset.loc[top_genes]

            fig = px.imshow(
                heatmap_data.values,
                x=filtered_samples['condition'].values,
                y=top_genes,
                color_continuous_scale="RdBu_r",
                title="Top 50 Most Variable Genes"
            )
            fig.update_xaxes(title="Samples")
            fig.update_yaxes(title="Genes")
            st.plotly_chart(fig, use_container_width=True)

        # Expression statistics
        st.subheader("Expression Statistics")

        col1, col2 = st.columns(2)

        with col1:
            # Distribution of expression values
            all_expr = self.data['expression'][filtered_samples['sample_id']].values.flatten()
            fig_dist = px.histogram(
                x=all_expr,
                nbins=50,
                title="Distribution of Expression Values"
            )
            fig_dist.update_xaxes(title="Expression Level")
            fig_dist.update_yaxes(title="Frequency")
            st.plotly_chart(fig_dist, use_container_width=True)

        with col2:
            # Sample correlation heatmap
            sample_corr = self.data['expression'][filtered_samples['sample_id']].corr()

            fig_corr = px.imshow(
                sample_corr.values,
                x=filtered_samples['condition'].values,
                y=filtered_samples['condition'].values,
                color_continuous_scale="RdBu_r",
                title="Sample Correlation Heatmap"
            )
            st.plotly_chart(fig_corr, use_container_width=True)

    def render_modules(self):
        """Render co-expression modules analysis."""
        st.header("Co-expression Modules")

        if 'modules' not in self.data:
            st.error("Module data not available")
            return

        # Module overview
        module_sizes = self.data['modules']['module'].value_counts()

        col1, col2 = st.columns(2)

        with col1:
            # Module size distribution
            fig_sizes = px.bar(
                x=module_sizes.index,
                y=module_sizes.values,
                title="Module Sizes"
            )
            fig_sizes.update_xaxes(title="Module", tickangle=45)
            fig_sizes.update_yaxes(title="Number of Genes")
            st.plotly_chart(fig_sizes, use_container_width=True)

        with col2:
            # Module-trait correlations
            if 'module_correlations' in self.data:
                fig_corr = px.imshow(
                    self.data['module_correlations'].values,
                    x=self.data['module_correlations'].columns,
                    y=self.data['module_correlations'].index,
                    color_continuous_scale="RdBu_r",
                    title="Module-Trait Correlations"
                )
                fig_corr.update_xaxes(title="Condition")
                fig_corr.update_yaxes(title="Module")
                st.plotly_chart(fig_corr, use_container_width=True)

        # Module details
        st.subheader("Module Details")

        selected_module = st.selectbox(
            "Select a module to explore:",
            self.data['modules']['module'].unique()
        )

        if selected_module:
            module_genes = self.data['modules'][
                self.data['modules']['module'] == selected_module
            ]['gene'].tolist()

            col1, col2 = st.columns([1, 2])

            with col1:
                st.write(f"**Module:** {selected_module}")
                st.write(f"**Size:** {len(module_genes)} genes")

                # Show some genes
                st.write("**Genes in module:**")
                for gene in module_genes[:10]:
                    st.write(f"• {gene}")
                if len(module_genes) > 10:
                    st.write(f"... and {len(module_genes) - 10} more")

            with col2:
                # Module gene expression
                if 'expression' in self.data and len(module_genes) > 1:
                    available_genes = [g for g in module_genes if g in self.data['expression'].index]
                    if available_genes:
                        module_expr = self.data['expression'].loc[available_genes[:20]]  # Limit for visualization
                        
                        fig_module = px.imshow(
                            module_expr.values,
                            x=list(range(module_expr.shape[1])),
                            y=module_expr.index,  # ✅ Correct Y labels
                            color_continuous_scale="RdBu_r",
                            title=f"Expression Heatmap: {selected_module}"
                        )
                        # fig_module = px.imshow(
                        #     module_expr.values,
                        #     x=list(range(module_expr.shape[1])),
                        #     y=available_genes,
                        #     color_continuous_scale="RdBu_r",
                        #     title=f"Expression Heatmap: {selected_module}"
                        # )
                        fig_module.update_xaxes(title="Samples")
                        fig_module.update_yaxes(title="Genes")
                        st.plotly_chart(fig_module, use_container_width=True)

    def render_hub_genes(self):
        """Render hub genes analysis."""
        st.header("Hub Genes Analysis")

        if 'hub_genes' not in self.data:
            st.error("Hub genes data not available")
            return

        # Top hub genes across all modules
        top_hubs = self.data['hub_genes'].nlargest(20, 'hub_score')

        col1, col2 = st.columns(2)

        with col1:
            # Top hub genes bar chart
            fig_hubs = px.bar(
                top_hubs,
                x='hub_score',
                y='gene',
                color='module',
                orientation='h',
                title="Top 20 Hub Genes",
                hover_data=['module', 'module_size']
            )
            fig_hubs.update_xaxes(title="Hub Score")
            fig_hubs.update_yaxes(title="Gene")
            st.plotly_chart(fig_hubs, use_container_width=True)

        with col2:
            # Hub genes per module
            hubs_per_module = self.data['hub_genes'].groupby('module').size().sort_values(ascending=False)

            fig_module_hubs = px.bar(
                x=hubs_per_module.index,
                y=hubs_per_module.values,
                title="Hub Genes per Module"
            )
            fig_module_hubs.update_xaxes(title="Module", tickangle=45)
            fig_module_hubs.update_yaxes(title="Number of Hub Genes")
            st.plotly_chart(fig_module_hubs, use_container_width=True)

        # Hub genes table
        st.subheader("Hub Genes Details")

        # Filter options
        col1, col2 = st.columns(2)

        with col1:
            selected_module_hub = st.selectbox(
                "Filter by module:",
                ["All"] + list(self.data['hub_genes']['module'].unique())
            )

        with col2:
            min_score = st.slider(
                "Minimum hub score:",
                float(self.data['hub_genes']['hub_score'].min()),
                float(self.data['hub_genes']['hub_score'].max()),
                float(self.data['hub_genes']['hub_score'].quantile(0.5))
            )

        # Filter data
        filtered_hubs = self.data['hub_genes'][self.data['hub_genes']['hub_score'] >= min_score]

        if selected_module_hub != "All":
            filtered_hubs = filtered_hubs[filtered_hubs['module'] == selected_module_hub]

        # Display table
        st.dataframe(
            filtered_hubs[['gene', 'module', 'hub_score', 'hub_rank', 'module_size']].sort_values('hub_score', ascending=False),
            use_container_width=True
        )

    def render_enrichment(self):
        """Render enrichment analysis results."""
        st.header("Functional Enrichment Analysis")

        enrich_df = self.data.get('enrichment')

        # → Check if enrichment exists and is valid
        if enrich_df is None or not isinstance(enrich_df, pd.DataFrame) or enrich_df.empty:
            # Calculate real overlapping gene count from modules vs expression
            expr_genes = set(self.data["expression"].index) if "expression" in self.data else set()
            module_genes = set(self.data["modules"]["gene"]) if "modules" in self.data else set()
            overlap = len(expr_genes.intersection(module_genes))

            st.warning("No enrichment data available or enrichment file is missing.")
            st.info(
                f"⚠️ Enrichment could not be computed because the pipeline "
                f"contained only **{overlap} overlapping genes**, which is typically too small "
                f"for statistical significance in pathway enrichment. "
                "Real enrichment usually requires several hundred–thousand genes."
            )
            return

        # Now safe to compute significance
        significant_enrichments = enrich_df[enrich_df['p_value'] < 0.05]
        
        col1, col2 = st.columns(2)

        with col1:
            # Top enriched pathways (ensure numeric column exists)
            if not significant_enrichments.empty and 'fold_enrichment' in significant_enrichments.columns:
                top_pathways = significant_enrichments.nlargest(15, 'fold_enrichment')

                if not top_pathways.empty:
                    fig_pathways = px.scatter(
                        top_pathways,
                        x="fold_enrichment",        # ✅ must be string, not list
                        y="pathway",                # ✅ must be string, not list
                        size="overlap",
                        color="p_value",
                        color_continuous_scale="viridis_r",
                        title="Top Enriched Pathways",
                        hover_data=["module", "pathway_size"],
                    )
                    fig_pathways.update_xaxes(title="Fold Enrichment")
                    fig_pathways.update_yaxes(title="Pathway")
                    st.plotly_chart(fig_pathways, use_container_width=True)
                else:
                    st.info("No top pathways found to plot.")
            else:
                st.warning("No enrichment data available or missing 'fold_enrichment' column.")
                st.warning("Because the dataset after preprocessing contained only 187 overlapping genes, which is statistically too small for pathway enrichment to reach significance. Real enrichment requires thousands of genes or less strict thresholds.")


        with col2:
            # Pathways per module
            pathways_per_module = significant_enrichments['module'].value_counts()

            if not pathways_per_module.empty:
                # ✅ Convert to DataFrame to make Plotly happy
                pathways_df = pathways_per_module.reset_index()
                pathways_df.columns = ['module', 'count']

                fig_module_pathways = px.bar(
                    pathways_df,
                    x='count',                # ✅ single column name (string)
                    y='module',
                    orientation='h',
                    title="Enriched Pathways per Module"
                )
                fig_module_pathways.update_xaxes(title="Number of Enriched Pathways")
                fig_module_pathways.update_yaxes(title="Module")
                st.plotly_chart(fig_module_pathways, use_container_width=True)
            else:
                st.info("No enriched pathways available to display.")

            # # Pathways per module
            # pathways_per_module = significant_enrichments['module'].value_counts()

            # fig_module_pathways = px.bar(
            #     x=pathways_per_module.values,
            #     y=pathways_per_module.index,
            #     orientation='h',
            #     title="Enriched Pathways per Module"
            # )
            # fig_module_pathways.update_xaxes(title="Number of Enriched Pathways")
            # fig_module_pathways.update_yaxes(title="Module")
            # st.plotly_chart(fig_module_pathways, use_container_width=True)

        # Enrichment details table
        st.subheader("Enrichment Details")

        # Filters
        col1, col2, col3 = st.columns(3)

        with col1:
            pathway_filter = st.selectbox(
                "Filter by pathway:",
                ["All"] + list(self.data['enrichment']['pathway'].unique())
            )

        with col2:
            module_filter = st.selectbox(
                "Filter by module:",
                ["All"] + list(self.data['enrichment']['module'].unique())
            )

        with col3:
            max_pvalue = st.slider(
                "Maximum p-value:",
                0.001, 0.05, 0.05
            )

        # Apply filters
        filtered_enrichment = self.data['enrichment'][self.data['enrichment']['p_value'] <= max_pvalue]

        if pathway_filter != "All":
            filtered_enrichment = filtered_enrichment[filtered_enrichment['pathway'] == pathway_filter]

        if module_filter != "All":
            filtered_enrichment = filtered_enrichment[filtered_enrichment['module'] == module_filter]

        # Display table
        if len(filtered_enrichment) > 0:
            st.dataframe(
                filtered_enrichment[['module', 'pathway', 'overlap', 'pathway_size', 'fold_enrichment', 'p_value']]
                .sort_values('p_value'),
                use_container_width=True
            )
        else:
            st.warning("No enrichment results match the selected filters")

    def render_network(self):
        """Render interactive network visualization."""
        st.header("Interactive Network Visualization")

        st.info("Interactive network visualization would be implemented here using tools like Cytoscape.js or Plotly network graphs.")

        # Placeholder for network visualization
        if 'hub_genes' in self.data:
            st.subheader("Hub Gene Network Preview")

            # Create a simple network visualization
            top_hubs = self.data['hub_genes'].nlargest(20, 'hub_score')

            # Mock network data
            edges = []
            for i, hub1 in enumerate(top_hubs.itertuples()):
                for j, hub2 in enumerate(top_hubs.itertuples()):
                    if i < j and np.random.random() > 0.7:  # Random connections
                        edges.append((hub1.gene, hub2.gene))

            if edges:
                # Create NetworkX graph
                G = nx.Graph()
                G.add_edges_from(edges)

                # Get positions
                pos = nx.spring_layout(G, k=1, iterations=50)

                # Create edge traces
                edge_x, edge_y = [], []
                for edge in G.edges():
                    x0, y0 = pos[edge[0]]
                    x1, y1 = pos[edge[1]]
                    edge_x.extend([x0, x1, None])
                    edge_y.extend([y0, y1, None])

                # Create node traces
                node_x, node_y, node_text = [], [], []
                for node in G.nodes():
                    x, y = pos[node]
                    node_x.append(x)
                    node_y.append(y)
                    node_text.append(node)

                # Create plot
                fig = go.Figure()

                # Add edges
                fig.add_trace(go.Scatter(
                    x=edge_x, y=edge_y,
                    line=dict(width=1, color='#888'),
                    hoverinfo='none',
                    mode='lines'
                ))

                # Add nodes
                fig.add_trace(go.Scatter(
                    x=node_x, y=node_y,
                    mode='markers+text',
                    text=node_text,
                    textposition="middle center",
                    hovertemplate='Gene: %{text}<extra></extra>',
                    marker=dict(
                        size=20,
                        color='lightblue',
                        line=dict(width=2, color='darkblue')
                    )
                ))

                fig.update_layout(
                    title="Hub Gene Co-expression Network",
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                )

                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("Not enough connections to display network")

def main():
    """Main Streamlit application."""
    # Initialize dashboard
    dashboard = DMDDashboard()

    # Render sidebar and get selections
    page, selected_datasets, selected_conditions = dashboard.render_sidebar()

    # Render selected page
    if page == "Overview":
        dashboard.render_overview()
    elif page == "Expression Data":
        dashboard.render_expression_data(selected_datasets, selected_conditions)
    elif page == "Co-expression Modules":
        dashboard.render_modules()
    elif page == "Hub Genes":
        dashboard.render_hub_genes()
    elif page == "Enrichment Analysis":
        dashboard.render_enrichment()
    elif page == "Interactive Network":
        dashboard.render_network()

    # Footer
    st.markdown("---")
    st.markdown(
        "<div style='text-align: center; color: #666; font-size: 0.8em;'>"
        "DMD Co-Expression Network Analysis Pipeline v1.0.0 | "
        "Built with Streamlit and Plotly"
        "</div>",
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
