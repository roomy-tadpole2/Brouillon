```python
import json
import pandas as pd

GB1_DATA = pd.read_csv(r'../data/SEC.csv', index_col=0, encoding='gbk')

INDEX = GB1_DATA.index
COLUMNS = GB1_DATA.columns

display(INDEX)
display(COLUMNS)
```


    Index(['Accuracy', 'macro_Precision', 'macro_Recall', 'macro_F1-score',
           'macro_AP', 'weighted_Precision', 'weighted_Recall',
           'weighted_F1-score', 'weighted_AP'],
          dtype='object', name='Metric')



    Index(['Pre-train-all', 'Fine-tune-all', 'Pre-train-prat', 'Fine-tune-part'], dtype='object')



```python
methods = []

for col in COLUMNS:
    if (col not in methods):
        methods.append(col)

values = {}
errors = {}

for method in methods:
    vals, errs = [], []
    for entry in GB1_DATA[method]:
        val, err = entry.split("±")
        vals.append(float(val))
        errs.append(float(err))
    values[method] = vals
    errors[method] = errs

display(methods)
display(values)
display(errors)
```


    ['Pre-train-all', 'Fine-tune-all', 'Pre-train-prat', 'Fine-tune-part']



    {'Pre-train-all': [0.829,
      0.822,
      0.834,
      0.827,
      0.901,
      0.831,
      0.829,
      0.83,
      0.904],
     'Fine-tune-all': [0.817,
      0.81,
      0.82,
      0.815,
      0.887,
      0.818,
      0.817,
      0.817,
      0.891],
     'Pre-train-prat': [0.845,
      0.835,
      0.844,
      0.84,
      0.909,
      0.845,
      0.845,
      0.845,
      0.913],
     'Fine-tune-part': [0.817,
      0.808,
      0.813,
      0.81,
      0.876,
      0.817,
      0.817,
      0.817,
      0.882]}



    {'Pre-train-all': [0.001, 0.001, 0.0, 0.001, 0.0, 0.0, 0.001, 0.001, 0.0],
     'Fine-tune-all': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     'Pre-train-prat': [0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     'Fine-tune-part': [0.001, 0.002, 0.001, 0.002, 0.0, 0.001, 0.001, 0.001, 0.0]}



```python
# pallete_name = 'pallete_new1'
pallete_name = 'pallete_new3'
pallete = json.load(open(f'../graphql/{pallete_name}.json', 'r'))
pallete = pallete[:len(methods)]
display(pallete)
```


    ['#ED7823', '#CCA763', '#D3D2D2', '#3B51A2']



```python
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(len(INDEX))
width = 0.16

fig, ax = plt.subplots(
    2, 1,
    figsize=(10, 4),
    gridspec_kw={
        'height_ratios': [5, 1],
        'hspace': 0.11
        }
    )

baseline = 0.75
buns = (0, 0.7)

for i, method in enumerate(methods):
    heights = np.array(values[method])
    heights -= baseline
    ax[0].bar(
        x + i * width,
        heights,
        bottom = baseline,
        width = width,
        yerr = errors[method],
        capsize = 2.5,
        label = method,
        color = pallete[i],
        edgecolor = "black",
        linewidth = 0.5,
        alpha = 0.75,
        error_kw = dict(ecolor="#444444", lw=0.5, capthick=0.5)
    )

for i, method in enumerate(methods):
    heights = [buns[1]] * len(INDEX)
    bottoms = [buns[0]] * len(INDEX)
    ax[1].bar(
        x + i * width,
        heights,
        bottom = bottoms,
        width = width,
        label = method,
        color = pallete[i],
        edgecolor = "black",
        linewidth = 0.5,
        alpha = 0.35,
    )

ax[0].set_ylabel("Score", fontsize=12)
ax[0].set_title("SEC Performance Comparison", fontsize=14)

ax[0].set_xticks(x + width * (len(methods) - 1) / 2)
ax[0].set_xticklabels([])

ax[0].set_ylim(baseline - 0.01, 0.95)

ax[0].spines["top"].set_visible(False)
ax[0].spines["right"].set_visible(False)
ax[0].spines["bottom"].set_visible(False)

ax[1].spines["top"].set_visible(False)
ax[1].spines["right"].set_visible(False)

#hide the xticks line in ax[0]
ax[0].tick_params(axis='x', which='both', length=0)

# complete copy the xticks and xlim from ax[0] to ax[1]
ax[1].set_xlim(ax[0].get_xlim())
ax[1].set_xticks(x + width * (len(methods) - 1) / 2)
ax[1].set_xticklabels(INDEX, fontsize=14, rotation=45, ha='right')

ax[1].set_ylim(buns[0], buns[1] + 0.1)

ax[1].set_yticks([0])

ax[0].legend(
    bbox_to_anchor=(1.05, 0.2),
    loc="center left",
    frameon=False,
    fontsize=14
)

for i in range(len(methods) - 1):
    for xpos in x + width * (i + 0.5):
        ax[0].text(
            xpos,
            -0.0,
            "⋮",
            transform=ax[0].get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=14,
            color="black"
        )

plt.savefig(rf"../figure/SEC/test_{pallete_name}.png", dpi=300, bbox_inches='tight')
plt.savefig(rf"../figure/SEC/test_{pallete_name}.pdf", bbox_inches='tight')
plt.show()
```


    
![png](plot_SEC_files/plot_SEC_3_0.png)
    

