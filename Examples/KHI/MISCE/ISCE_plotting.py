import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# These are defaults for screen plotting, not paper
full_width = (7.68, 4.8)
half_width = (3.84, 4.8)
half_double = (7.68, 9.6)
plt.rcParams['figure.figsize'] = half_width # 50% of the screen, vertical, assuming a 1x2 plot
plt.rcParams['font.size'] = 6
save_plots = False

var_dict = {
    'n': { 'group': 'Primitive', 'legend': r"$n$"},
    'rho': { 'group': 'Primitive', 'legend': r"$\rho$"},
    'v1': { 'group': 'Primitive', 'legend': r"$v_x$"},
    'v2': { 'group': 'Primitive', 'legend': r"$v_y$"},
    'v3': { 'group': 'Primitive', 'legend': r"$v_z$"},
    'p': { 'group': 'Primitive', 'legend': r"$p$"},
    'Pi': { 'group': 'Primitive', 'legend': r"$\Pi$"},
    'q1': { 'group': 'Primitive', 'legend': r"$q_{x}$"},
    'q2': { 'group': 'Primitive', 'legend': r"$q_{y}$"},
    'q3': { 'group': 'Primitive', 'legend': r"$q_{z}$"},
    'pi11': { 'group': 'Primitive', 'legend': r"$\pi_{xx}$"},
    'pi12': { 'group': 'Primitive', 'legend': r"$\pi_{xy}$"},
    'pi13': { 'group': 'Primitive', 'legend': r"$\pi_{xz}$"},
    'pi22': { 'group': 'Primitive', 'legend': r"$\pi_{yy}$"},
    'pi23': { 'group': 'Primitive', 'legend': r"$\pi_{yz}$"},
    'pi33': { 'group': 'Primitive', 'legend': r"$\pi_{zz}$"},
    'T': {'group': 'Auxiliary', 'legend': r"$T$"},
    'W': {'group': 'Auxiliary', 'legend': r"$W$"},
    'h': {'group': 'Auxiliary', 'legend': r"$h$"},
    'PiNS': {'group': 'Auxiliary', 'legend': r"$\Pi_{NS}$"},
    'a1': { 'group': 'Auxiliary', 'legend': r"$a_{x}$"},
    'a2': { 'group': 'Auxiliary', 'legend': r"$a_{y}$"},
    'a3': { 'group': 'Auxiliary', 'legend': r"$a_{z}$"},
}

def plot_khi(fss, group_str = 'Primitive', quant_str='n', quant_legend_str = r"$n$", 
             model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$\eta$", figsize=full_width):
    
    assert len(fss) % 2 == 0

    access_str = group_str + '/' + quant_str

    coeff = fss[0]['Optional'].attrs[coeff_str][0]
    xmin = fss[0]['Domain'].attrs['xmin'][0]
    xmax = fss[0]['Domain'].attrs['xmax'][0]
    ymin = fss[0]['Domain'].attrs['ymin'][0]
    ymax = fss[0]['Domain'].attrs['ymax'][0]
    extent = [ymin, ymax, xmin, xmax]
    max_quant = -np.inf
    min_quant = np.inf
    for fs in fss:
        max_quant = max(max_quant, fs[access_str][:].max())
        min_quant = min(min_quant, fs[access_str][:].min())

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(len(fss)//2, 4, width_ratios=[9, 9, 1.5, 1], wspace=0.01, hspace=0.05)
    for i, fs in enumerate(fss):
        ax = fig.add_subplot(gs[i//2, i%2])
        im = ax.imshow(fs[access_str][:], vmin=min_quant, vmax=max_quant, extent=extent)
        #if i == 0:
            #ax.set_title(rf"{model_str}, {quant_legend_str}")
        #if i == 1:
            #ax.set_title(rf"{coeff_legend_str}={coeff:.2g}")
        ax.set_ylabel(rf"$t = ${fs.attrs['t'][0]:.2f}")
        if i%2 == 1:
            ax.yaxis.set_label_position('right')
        ax.set_xticks([])
        ax.set_yticks([])
    ax_cb1 = fig.add_subplot(gs[:, -1])
    fig.colorbar(im, cax=ax_cb1)

    return fig

def plot_shock(fss, plot_shape=None, group_str='Primitive', quant_str='n', quant_legend_str = r"$n$", 
               model_str='MISCE', coeff_str='eta', coeff_legend_str = r"$\eta$", figsize=full_width):

    access_str = group_str + '/' + quant_str

    coeff = fss[0]['Optional'].attrs[coeff_str][0]
    xmin = fss[0]['Domain'].attrs['xmin'][0]
    xmax = fss[0]['Domain'].attrs['xmax'][0]
    x = np.linspace(xmin, xmax, fss[0]['Domain'].attrs['nx'][0])
    max_quant = -np.inf
    min_quant = np.inf
    for fs in fss:
        max_quant = max(max_quant, fs[access_str][:].max())
        min_quant = min(min_quant, fs[access_str][:].min())
    dq = max_quant - min_quant
    if dq < 1e-6:
        if abs(max_quant) < 1e-6:
            dq = 1
        else:
            dq = abs(max_quant)
    q_max = max_quant + 5e-2*dq
    q_min = min_quant - 5e-2*dq

    if plot_shape is not None:
        fig, axes = plt.subplots(plot_shape[0], plot_shape[1], sharex=True, sharey=True,
                                 figsize=figsize, gridspec_kw={'wspace':0.01, 'hspace':0.001})
    else:
        fig, axes = plt.subplots(1, len(fss), sharex=True, sharey=True,
                                 figsize=figsize, gridspec_kw={'wspace':0.01, 'hspace':0.001})
    
    for i, (fs, ax) in enumerate(zip(fss, axes.flatten())):
        ax.plot(x, fs[access_str][:].flatten())
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(q_min, q_max)
    return fig

def plot_shock_compare(fss, 
                       group_strs=('Primitive',), quant_strs=('n',), quant_legend_strs = (r"$n$",), 
                       model_strs=('Ideal', 'MISCE',), coeff_strs=(None,'zeta'), 
                       coeff_legend_strs = (None, r"$\zeta$"), 
                       figsize=full_width, lss=('-', '--')):

    assert len(fss) == 2
    assert len(model_strs) == 2
    assert len(quant_strs) % 2 == 0
    fig, axes = plt.subplots(len(quant_strs)//2, 2,
                             figsize=figsize, gridspec_kw={'wspace':0.01, 'hspace':0.001})
    
    for i, ax in enumerate(axes.flatten()):
        access_str = group_strs[i]+'/'+quant_strs[i]
        max_quant = -np.inf
        min_quant = np.inf
        for model in range(2):
            xmin = fss[model]['Domain'].attrs['xmin'][0]
            xmax = fss[model]['Domain'].attrs['xmax'][0]
            x = np.linspace(xmin, xmax, fss[model]['Domain'].attrs['nx'][0])
            max_quant = max(max_quant, fss[model][access_str][:].max())
            min_quant = min(min_quant, fss[model][access_str][:].min())

            label_str = model_strs[model]
            if coeff_strs[model] is not None:
                coeff = fss[model]['Optional'].attrs[coeff_strs[model]][0]
                label_str += ", "+coeff_legend_strs[model]+f"={coeff:.2g}"
            ax.plot(x, fss[model][access_str][:].flatten(), ls=lss[model], label=label_str)
        ax.set_xlim(xmin, xmax)
        dq = max_quant - min_quant
        if dq < 1e-6:
            if abs(max_quant) < 1e-6:
                dq = 1
            else:
                dq = abs(max_quant)
        q_max = max_quant + 5e-2*dq
        q_min = min_quant - 5e-2*dq
        ax.set_ylim(q_min, q_max)
        ax.set_ylabel(quant_legend_strs[i])
        if i == 1:
            ax.legend()
        if i%2 == 1:
            ax.yaxis.set_label_position('right')
            ax.yaxis.tick_right()
        if not (i == len(quant_strs)-2):
            ax.set_xticks([])
    return fig

def plot_shock_compare_wrapper(fss, 
                       quant_strs=('n',),
                       model_strs=('Ideal', 'MISCE',), coeff_strs=(None,'zeta'),
                       figsize=full_width, lss=('-', '--')):
    coeff_legend_strs = []
    for coeff_str in coeff_strs:
        if coeff_str is not None:
            coeff_legend_strs.append(rf"$\{coeff_str}$")
        else:
            coeff_legend_strs.append(None)
    group_strs = [var_dict[var]['group'] for var in quant_strs]
    quant_legend_strs = [var_dict[var]['legend'] for var in quant_strs]
    return plot_shock_compare(fss, group_strs, quant_strs, quant_legend_strs, model_strs, 
                              coeff_strs, coeff_legend_strs, figsize, lss)

# fig=plot_shock_compare(scfs, group_strs=('Primitive', 'Primitive', 'Auxiliary', 'Auxiliary'), quant_strs=('v1', 'rho', 'PiNS', 'T'), model_strs=('Ideal', 'Bulk'), quant_legend_strs=(r"$v_x$", r"$\rho$", r"$\Pi_{NS}$", r"$T$"))
