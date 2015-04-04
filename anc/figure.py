# matplotlib parameters
fparams = {'axes.labelsize': 22,
           'axes.titlesize': 20,
           'text.fontsize': 18,
           'legend.fontsize': 14,
           'axes.linewidth': 1.5,
           'font.family': 'serif',
           'font.serif': 'Computer Modern Roman',
           'xtick.labelsize': 20,
           'xtick.major.size': 5.5,
           'xtick.major.width': 1.5,
           'ytick.labelsize': 20,
           'ytick.major.size': 5.5,
           'ytick.major.width': 1.5,
           'text.usetex': True,
           'figure.autolayout': True}
plt.rcParams.update(fparams)


f, ax = plt.subplots(figsize=(7, 4))
#(6, 5) (7.7, 4.2) (7.5, 4) (6, 4) (4, 3) (9, 4.7)

im = ax.imshow(dos, vmin=0, vmax=150,
               interpolation='none',
               extent=(0, 2*np.pi, 0, 1),
               aspect='auto',
               cmap='gist_heat_r')

ax.plot(a, b, c=col, marker='.', markersize=4.5, linewidth=0.5)
ax.plot(a, b, linewidth=2, ls='--', color='black')
ax.scatter(a, b, marker='.', color='k')

ax.set_title(r'', y=1.08, x=0.58, fontsize=16)
ax.yaxis.set_label_position('right')
ax.yaxis.tick_right()

ax.set_xlabel(r'', labelpad=-8, x=0.25, fontsize=24)
ax.set_ylabel(r'', labelpad=-5, y=0.25, rotation=0)

ax.set_xticks([0., np.pi, 2*np.pi])
ax.set_xticklabels(['$0$', '$\pi$', '$2\pi$'])
#ax.set_xticklabels([])
ax.set_xlim([0, 2*np.pi])
#same with x->y

cbar = f.colorbar(im, cax=cbaxes, orientation='horizontal')
cbar.set_ticks([0, 1])
cbar.set_ticklabels(['$0$', '$\Delta$'])
cbar.set_label('$\epsilon$', rotation=0, labelpad=-15, y=0.5)
cbar.solids.set_edgecolor("face")

f.savefig('f.pdf', transparent=True, bbox_inches='tight', pad_inches=0.05)
plt.close(f)


#ax.set_aspect(3)


f, axes = plt.subplots(3, 1, figsize=(2.5, 7.3))
f.subplots_adjust(hspace=0.17)

for i, (ax, level, title) in enumerate(zip(axes, levels, titles)):
    im = ax.imshow(..)

    if i == 1:
        ax.plot(..)
    elif i == 2:
        ax.plot(..)

f.subplots_adjust(top=0.9)
cbaxes = f.add_axes([0.2, 0.02, 0.6, 0.015])
