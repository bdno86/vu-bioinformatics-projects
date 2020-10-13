# ROC Assigment

## Question 2 - 10 pts
To evaluate a diagnostic method, the probability that the test marks a healthy person as sick is evaluated. This is called the False Positive Rate (FPR). Note that this figure does not depend on the a priori distribution of sick and healthy persons. 

$FPR\:=\:\frac{fp}{fp+tn}$

Where fp is all the false positives and tn all the true negatives.

This probability of "false alarms" is not the only thing to be worried about. The probability that a test marks a healthy person as healthy is also important. This is called the True Positive Rate (TPR). It can be expressed as a formula as:

$TPR\:=\:\frac{tp}{tp+fn}$

Here the tp is all the true positives and fn all the false negatives. 

Because it is easy to get a low FPR and easy to get a low TPR, but not easy to do both, it is clear that some trade-off exist. To find out how well a method separates the positives from the negatives, a ROC-plot is a useful tool. A parameter of a classifier method is the discrimination threshold, which determines at which level the entities to be classified are called similar or different. Varying this threshold yields different TPR and FPR predicted by the classifier. To create a ROC plot, the discrimination threshold is varied and for each threshold value, the FPR is plotted against the TPR rate. 

The performance of tools that find homologous of proteins can be assessed in the same way. We will define a protein pair that is homologous as a positive, and a protein pair that is different as a negative. If a protein-pair is correctly predicted to be homologous by BLAST, it is called a true positive. If BLAST predicts a protein-pair as homologous, while the proteins are not related, it is called a false positive. Of course it is not always clear if a protein is homologous.

**A) There is a trivial method that yields a FPR of zero. Please state such a method below. [5]**

A trivial way to achieve an FPR of zero is by selecting a classifier threshold so that none of the instances will be classified as positive at all. When everything is classified as negative, nothing can be either true or false positive.

Given:
$$FPR=\frac{FP}{FP+TN}$$ 
Then:
$$given \quad TN>0 \quad and \quad FP=0:\\ \space \\
FPR=\frac{FP}{FP+TN}=\frac{0}{0+TN}=\frac{0}{TN}=0
$$

**B) There is a trivial method that yields a TPR of one. Please describe it. [5]**

A trivial way to achieve an TPR of one is by selecting a classifier threshold so that all of the instances will be classified as positive. When everything is classified as positive, nothing will be classified as negative and therefor there can't be any true or false negatives.

Given:
$$TPR=\frac{TP}{TP+FN}$$ 
Then:
$$given \quad TP>1 \quad and \quad FN=0: \\ \space \\
TPR=\frac{TP}{TP+FN}=\frac{TP}{TP+0}=\frac{TP}{TP}=1
$$

## Question 3 - 10 pts
Methods usually provide a prediction in the form of a numeric value. We can use this value together with a varying threshold, to assign positive and negative predictions. In this way we can calculate the TPR and the FPR. A resulting ROC plot visualises the trade-off that exists in minimizing the FPR while maximizing the TPR.

**A) What would the ROC-plot of a method that works randomly look like? [5]**

A method or model that classifies any instance given to it by random will generate a plot that tends toward a diagonal line. With a large enough sample size the line will be a completely straight diagonal line.

The AUC (Area Under Curve) is equal to 0.5.

An example of a perfectly random method from which a very large number of samples is drawn is shown in the graphic below indicated by C.


**B) What would the ROC-plot of a method that perfectly separates the two groups look like? [5]**

A perfect predictor will produce a perfect TPR of 1 and will have no false positives at all, yielding an FPR of 0. For this reason, the plot will be a straight vertical line from (0.0, 0.0) to (0, 1.0) followed by a horizontal line from (0,1.0) to (1.1, 1.0). The AUC is equal to 1.0.

Or simply put: straight up and then all the way to the right.

An example of a perfect method from which a very large number of samples is drawn is shown in the graphic below indicated by A.

![enter image description here](https://www.ahajournals.org/cms/asset/4a1ab89e-38b2-4c50-a4fb-5f64be9f868d/16ff2.jpeg)
https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.105.594929
Zou, K. H., O’Malley, A. J., & Mauri, L. (2007). Receiver-operating characteristic analysis for evaluating diagnostic tests and predictive models. Circulation, 115(5), 654-657.


## Question 4 10 pts
Finish the skeleton scriptPreview the document to create ROC curves. The ROC curves should indicate the performance of BLAST on the GO database. Use your results from the last two practicals. 

```sh
python3 roc_plot_skeleton.py -blast_results BLAST_RESULTS.txt -go_results GO_RESULTS.txt -outpng 
FILENAME_ROC_PLOT.png
```
DONE

## Question 5 - 20 pts
**A) Show the generated ROC plot benchmarking BLAST against GO here. Please indicate all the parameters (e.g. e-value threshold, database used, GO score cut-off) used to generate this plot.**

**B) What can you conclude about the quality from your BLAST predictions?**

## Question 6 - 20 pts
**A) Also show a ROC plot for the PSI-BLAST runs, using similar setting as above. Please indicate all the parameters (e.g. e-value threshold, database used, GO score cut-off) used to generate this plot.**

**B) Do you see a difference in the ROC-curves between BLAST and PSI-BLAST? If yes, why and if no, why not? Please discuss.**

## Question 7 - 20 pts
**A) Which pair of proteins that is classified as different in GO has the lowest e-value?**

Q9UWN7,P13393 different 1,32e-45

Method we used to find the different pair with the lowest e-value:
```shell

join ./go_output_different_without_scores.txt output.txt | sort -g --key=3 // GREP EXCLUDE NA

```

**B) Where do you find this pair of proteins in the ROC-plot? What is the e-value?**

The e-value of this 
1,32e-45

**C) Do you think it is likely that these proteins are not homologous?**

**Explain your findings.**

## Question 8 - 10 pts
**Do you think the GO database is suitable as a gold standard to test the performance of PSI-BLAST? Keep in mind what PSI-BLAST is generally used for and what it actually does.**