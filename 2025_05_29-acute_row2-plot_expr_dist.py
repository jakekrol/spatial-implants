#!/usr/bin/env python3

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.stats import linregress

parser = argparse.ArgumentParser(description='Plot expression over distance for gene-(sets) of acute row2')
parser.add_argument('-i', '--input', type=str, required=True, help='Input table with expression and distance data')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file for the plot')
parser.add_argument('-f', '--fontsize', type=int, default=25, help='Font size for the plot')
args = parser.parse_args()

sizeannotate = args.fontsize - 5

df = pd.read_csv(args.input, sep='\t')
df.columns = ['gene', 'dist', 'expr', 'class']
# drop distance<=0
df = df[df['dist'] > 1]

# group by gene and plot each gene's expression over distance
df_grouped = df.groupby('gene')
for gene, group in df_grouped:
    # do linear regression
    slope, intercept, r_value, p_value, std_err = linregress(group['dist'], group['expr'], alternative='less')
    # plot points
    plt.scatter(group['dist'], group['expr'], label=gene)
    # plot regression line
    x_fit = pd.Series(range(2, group['dist'].max() + 1))
    y_fit = intercept + slope * x_fit
    # annotate p-value and r-squared
    plt.annotate(f'p={p_value:.2e}, r^2={r_value**2:.2f}', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=sizeannotate,
                 bbox=dict(facecolor='white', alpha=0.5))
    plt.plot(x_fit, y_fit, color='red', linestyle='--')
    plt.xlabel('Distance (spots away from implant)', fontsize=args.fontsize)
    plt.ylabel('Average expression', fontsize=args.fontsize)
    plt.title(gene, fontsize=args.fontsize)
    plt.tight_layout()
    plt.savefig(f"{args.output}_{gene}.png")
    plt.close()
    print(f"Plot saved for gene: {gene}")
