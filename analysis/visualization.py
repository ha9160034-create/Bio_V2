# 3D molecular visualization using py3Dmol and SASA comparative Plotly charts
import py3Dmol
import plotly.graph_objects as go
import pandas as pd

def render_protein_3d(pdb_string, bg_color='#111', style_type='cartoon',
                      show_surface=True, surface_opacity=0.3, surface_type='MS',
                      mutations=None, mut_color='red', zoom_to_mutations=False,
                      focus_mut=None, show_ligands=True,
                      ligand_style='Sticks & Spheres'):
    """Generate HTML for 3D protein structure visualization using py3Dmol."""
    view = py3Dmol.view(width="100%", height=450)
    view.addModel(pdb_string, 'pdb')
    view.setBackgroundColor(bg_color)

    # Set protein display style (cartoon, stick, sphere)
    style_dict = {style_type: {'color': 'spectrum'}}
    view.setStyle({'model': -1}, style_dict)

    # Display ligand and heteroatom (HETATM) molecules
    if show_ligands and ligand_style != "Hide":
        if "Sticks Only" in ligand_style or ligand_style == "sticks":
            het_style = {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.35}}
        elif "Spheres" in ligand_style and "Sticks" not in ligand_style:
            het_style = {'sphere': {'colorscheme': 'greenCarbon', 'radius': 1.0}}
        else:
            het_style = {
                'stick': {'colorscheme': 'greenCarbon', 'radius': 0.35},
                'sphere': {'colorscheme': 'greenCarbon', 'radius': 0.6}
            }
        view.addStyle({'hetflag': True}, het_style)
    else:
        view.setStyle({'hetflag': True}, {})

    # Display protein outer surface (applied only to protein residues, excluding ligands)
    if show_surface:
        surf_kind = py3Dmol.MS if surface_type == 'MS' else (py3Dmol.VDW if surface_type == 'VDW' else py3Dmol.SAS)
        view.addSurface(surf_kind, {'opacity': surface_opacity, 'color': '#FFC107'}, {'hetflag': False})

    # Color and highlight mutation sites (merge styles into one dict to prevent overwriting)
    if mutations:
        for mut in mutations:
            view.addStyle(mut, {
                style_type: {'color': mut_color},
                'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.3},
                'sphere': {'color': mut_color, 'radius': 1.2}
            })

    # Camera zoom control
    if focus_mut:
        view.zoomTo(focus_mut)
    elif zoom_to_mutations and mutations:
        view.zoomTo({'or': mutations})
    else:
        view.zoomTo()

    style_override = "<style>html, body { margin: 0 !important; padding: 0 !important; overflow: hidden !important; width: 100% !important; height: 100% !important; } div[id^='3dmolviewer'] { margin: 0 !important; padding: 0 !important; overflow: hidden !important; width: 100% !important; height: 100% !important; }</style>"
    return style_override + view._make_html()

def build_sasa_figure(comparison_df: pd.DataFrame) -> go.Figure:
    """Build an interactive Plotly chart comparing SASA (Solvent-Accessible Surface Area) profiles and highlighting mutations."""
    fig = go.Figure()

    # Plot Wild-Type reference curve
    fig.add_trace(go.Scatter(
        x=comparison_df['plot_index'],
        y=comparison_df['SASA_H'],
        name='Wild-Type (WT — Healthy Reference)',
        line=dict(color='#00ff88', width=2),
        fill=None,
        customdata=comparison_df[['res_num_h', 'WT_Res']].values,
        hovertemplate="🟢 WT (Wild-Type Reference): %{customdata[1]} (Residue %{customdata[0]}) — SASA (Solvent-Accessible Surface Area): %{y:.2f} Å² (square Angstroms)<extra></extra>"
    ))

    # Plot Mutant variant curve with shaded difference
    fig.add_trace(go.Scatter(
        x=comparison_df['plot_index'],
        y=comparison_df['SASA_M'],
        name='Mutant (MT — Variant Structure)',
        line=dict(color='#ff3333', width=2),
        fill='tonexty', 
        fillcolor='rgba(255, 51, 51, 0.1)',
        customdata=comparison_df[['res_num_m', 'MT_Res']].values,
        hovertemplate="🔴 MT (Mutant Variant): %{customdata[1]} (Residue %{customdata[0]}) — SASA (Solvent-Accessible Surface Area): %{y:.2f} Å² (square Angstroms)<extra></extra>"
    ))

    # Highlight mutation sites with yellow stars
    mutations_df = comparison_df[comparison_df['WT_Res'] != comparison_df['MT_Res']]
    if not mutations_df.empty:
        fig.add_trace(go.Scatter(
            x=mutations_df['plot_index'],
            y=mutations_df['SASA_M'],
            mode='markers',
            name='Mutation Sites',
            marker=dict(color='yellow', size=10, symbol='star', line=dict(color='black', width=1)),
            hovertemplate="⭐ Mutation: %{customdata[2]} (%{customdata[0]}) → %{customdata[3]} (%{customdata[1]})<br>Impact: %{customdata[4]}<extra></extra>",
            customdata=mutations_df[['res_num_h', 'res_num_m', 'WT_Res', 'MT_Res', 'Impact']].values
        ))

    # Use actual residue numbers as x-axis tick labels
    tick_step = max(1, len(comparison_df) // 20)
    fig.update_layout(
        template="plotly_dark",
        height=450,
        hovermode="x unified",
        xaxis=dict(
            title="Residue Position",
            rangeslider=dict(visible=True),
            tickmode='array',
            tickvals=comparison_df['plot_index'][::tick_step].tolist(),
            ticktext=comparison_df['res_num_label'][::tick_step].astype(str).tolist()
        ),
        yaxis=dict(title="SASA (Solvent-Accessible Surface Area), Å² (square Angstroms)"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=0, r=0, t=30, b=0)
    )
    return fig
