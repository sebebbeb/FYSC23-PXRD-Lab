import powerxrd as xrd
import matplotlib.pyplot as plt
import os
import re

# Monkey-patch the Chart class __init__ to override default K and lambdaKa
original_init = xrd.Chart.__init__

def new_init(self, x, y):
    '''
    Change the default K and lambdaKa values for the Chart class to fit 
    what the FYSC23 PXRD course uses.
    '''
    original_init(self, x, y)
    self.K = 0.94
    self.lambdaKa = 0.07107 # nm

xrd.Chart.__init__ = new_init

def backsub_multiplt():
    fig, axs = plt.subplots(2, 1, figsize=(6,8), sharex=True)
    fig.subplots_adjust(hspace=0.3)

    for ax in axs:
        ax.tick_params(labelbottom=True)

    # xrd.Data import tab-separated files (.xy) file with importfile()
    for i in range(2):
        data = xrd.Data('samples/filtered_sample{}.xy'.format(i+1)).importfile()
        chart = xrd.Chart(*data)
        axs[i].plot(*chart.backsub(), color='k', label='Sample {}'.format(i+1))
        axs[i].legend()
        axs[i].set_xlabel('2 $\\theta$ (deg)')
        axs[i].set_ylabel('Intensity (a.u.)')
        axs[i].grid()

    axs[0].set_title(f'PXRD Data With Background Subtraction for Sample 1 and 2')
    plt.savefig('./figs/backsub_multiplt.pdf')
    plt.show()

def all_peaks(filename):
    data = xrd.Data(filename).importfile()
    chart = xrd.Chart(*data)
    chart.backsub(tol=1.0, show=True)
    chart.allpeaks(tols=(0.1, 0.8), verbose=True, show=True)
    if filename == r'samples\filtered_sample1.xy':
        chart.SchPeak(xrange=[34.5,36.5], verbose=True, show=True)
    else:
        pass
    
    plt.grid()
    plt.xlabel('2 $\\theta$ (deg)')
    
    base = os.path.basename(filename)
    match = re.search(r'sample(\d+)', base, re.IGNORECASE)
    if match:
        sample_num = match.group(1)
    else:
        sample_num = os.path.splitext(base)[0]
    plt.suptitle(f'Backsub & Scherrer Width Calculation of All Peaks for Sample {sample_num}')
    plt.savefig(f'./figs/allpeaks_sample{sample_num}.pdf')
    plt.show()

def test_sch():
    data = xrd.Data('samples/filtered_sample1.xy').importfile()
    chart = xrd.Chart(*data)

    chart.backsub(tol=1.0, show=True)
    chart.SchPeak(xrange=[35, 36], verbose=True, show=True)
    plt.xlabel('2 $\\theta$')
    plt.title('backsub and Scherrer width calculation')
    plt.show()

if __name__ == '__main__':
    backsub_multiplt()
    all_peaks(r'samples\filtered_sample1.xy')
    # test_sch()
