# Assignment GO

## Notes

Go-terms for a single protein can be retrieved by accessing `http://www.uniprot.org/uniprot/?query=accession:<PROTID>&format=tab&columns=id,organism-id,go-id,pathway` and by replacing `<PROTID>` with the protein ID.

Command for running the script of this assignment:
```sh
python3 classify_go.py -uniprot SCOP_selections.txt -threshold1 0.111 -threshold2 0.054 -output_file results/go_output.txt
```

## Question 2
**Q:** As the name suggests, GO is an ontology. An ontology is an explicit specification of how terms in a specific domain are organised. in GO, the domain consists of terms annotating genes and gene products. These terms are organised in a specific way. This question consists out of two (A and B) questions.  

**A.** How are the terms in the Gene Ontology database organised? Be specific. (10 points)
**Answer:** 
Ontologies consist of objects or concepts that are connected to each other and describe their properties through relationships. These relationships form a graph that can be hierarchical and is often represented in graphical form.

The Gene Ontology describes gene products in terms of three aspects through three separate (but interrelated through some relation types) ontologies:
1. biological processes
2. molecular functions
3. cellular components

The goals of the GO term is describing the functions, processes and locations of the gene products. 
Every GO term has four properties:
1. Gene or gene product identifier
2. Go term ID
3. Reference code
4. Evidence code

The ontology structure is based on relationships. The following relationships can be distinguished:
1. Is_a
2. part_of
3. Regulates
	3a. postively_regulates
	3b. Negatively_regulates
4. Has_part
5. occurs_in

All together, the GO ontology can be visualized as an acyclic graph, where all the terms represent relationships to one another. 

**B.** Can you have multiple GO terms per protein? Are GO terms unique to only one protein? Explain why or why not. (10 points)
**Answer:** 
It is possible to have multiple GO terms per protein. GO terms can describe all kinds of properties and relations, and it makes sense that any single protein has a plethora of properties that can each individually be represented as a GO term.

GO terms are not unique to a single protein. GO terms often describe properties or subclasses of other information units, and therefore by their very nature cannot be exclusive to any single protein (unless you introduce an enormous amount of redundancy which would defeat the purpose of such a database).

## Question 3

**Q:** Biologists commonly use sequence alignment to infer whether two proteins are functionally similar. Typically, a sequence identity larger than 35% is considered to be a strong indication for functional similarity. On the other hand, a sequence identity between 20% and 35% corresponds to the so called twilight zone, for which there is not enough evidence to infer functional similarity. 
It might be the case that two proteins share exactly the same function, but have different sequences. Give a biological reason for why this may be the case. 
**Answer:**
A biological reason for this lies in the environment or the environmental factors multiple organisms (that are not ancestrally related to each other) are subjected to and subsequently have to adapt in similar ways to cope with those factors.

When in different settings same evolutionary pressures are exerted to different organisms, similar mechanisms and biological outcomes can be reached despite very different underlying genetic codes. The term for this is evolutionary convergence. This is also despite having no related common ancestors.

Evolutionary divergence, however, is a mechanism through which proteins/organisms start developing differently and genetically start to diverge from their ancestor. 

## Question 4
**Q:** An alternative way to infer functional similarity is by comparing GO terms. For example, proteins LAP3 (Uniprot ID: Q01532) and BLMH (Uniprot ID: Q13867) are functionally related. You may use QuickGo (Links to an external site.).

List the GO terms of these two proteins. Also list the GO terms that these two proteins have in common. What is the most specific similarity of these proteins?

**Answer:** 

*Go terms they have in common:*
1. GO:0005737 cytoplasm
2. GO:0008234 cysteine-type peptidase activity
3. GO:0006508 proteolysis
4. GO:0004197 cysteine-type endopeptidase activity
5. GO:0008233 peptidase activity
6. GO:0009636 response to toxic substance
7. GO:0016787 hydrolase activity

*Most specific similarity:*

As opposed to the answer in question 6 for which we used the uniprot database, here we used the GO terms as found in the quickGO database.

The most specific similarity between these two proteins is that they both have cysteine-type endopeptidase activity as indicated by their common designated GO term GO:0004197 - cysteine-type endopeptidase activity.

This GO term being the most specific one is easily deduced by looking at the other intersectional GO terms, which are either less specific forms of peptidase activity or very general activity (response to toxic substance).

## Question 5

**Q:** The idea is that the more GO terms two protein share, the more likely it is that they have similar functions. We can design a score taking this aspect into account. Preferably the score should be normalized, eg. between 0 and 1. 

Why is it convenient that in general scores are normalized? 

**Answer:** 
Normalization of scores is convenient because it allows values to be compared on a common scale. 

When a score is added to the group of scores, this one new score can be compared to any score already in the dataset and the distance between the two can be compared relative to the entire set on a common scale.

## Question 6

Q: Given a protein pair (p1, p2), a possible scoring function is: 

$$s(p1, p2) = \begin{Bmatrix}
\frac{\left | g(p1) \bigcap g(p2) \right |}{\left | g(p1) \bigcup g(p2) \right |} & \textrm{if } g(p1) \neq  \varnothing \textrm{ or } g(p2) \neq  \varnothing  \\ 
0 & \textrm{if } g(p1) = \varnothing \textrm{ and } g(p2) = \varnothing 
\end{Bmatrix}$$

where g(p) is the set of GO terms associated to protein p.  
(see also https://en.wikipedia.org/wiki/Union_(set_theory) (Links to an external site.) and https://en.wikipedia.org/wiki/Intersection_(set_theory) (Links to an external site.) )

**A)** Use the formula to compute the score between LAP3 and BLMH. (5 points)

**B)** What would be the score of a protein when compared to itself? (5 points) 

### PART A

Go-terms for protein **LAP3** - Q01532:
```
Entry	Organism ID	Gene ontology IDs	Pathway
Q01532	559292	GO:0000122; GO:0000978; GO:0003690; GO:0003697; GO:0003729; GO:0004197; GO:0005737; GO:0005739; GO:0008234; GO:0009636; GO:0043418; GO:0046677	
```

Go-terms for protein **BLMH** - Q13867:
```
Entry	Organism ID	Gene ontology IDs	Pathway
Q13867	9606	GO:0000209; GO:0004177; GO:0004180; GO:0004197; GO:0005634; GO:0005737; GO:0005829; GO:0006508; GO:0008234; GO:0009636; GO:0042493; GO:0042802; GO:0043418; GO:0070062	
```

lap3 $\cap$ blmh:
```python
>>> lap3 = set("GO:0000122; GO:0000978; GO:0003690; GO:0003697; GO:0003729; GO:0004197; GO:0005737; GO:0005739; GO:0008234; GO:0009636; GO:0043418; GO:0046677".split(';')) 
>>> blmh = set("GO:0000209; GO:0004177; GO:0004180; GO:0004197; GO:0005634; GO:0005737; GO:0005829; GO:0006508; GO:0008234; GO:0009636; GO:0042493; GO:0042802; GO:0043418; GO:0070062".split(';'))

# intersection = g(p1) intersection g(p2)
>>> intersection = lap3.intersection(blmh)
{' GO:0004197', ' GO:0005737', ' GO:0008234', ' GO:0009636', ' GO:0043418'}
# len(intersection)
>>> len(intersection)
5

# union = g(p1) union g(p2)
>>> union = lap3.union(blmh)
{' GO:0005737', ' GO:0005634', ' GO:0070062', ' GO:0004180', ' GO:0005829', ' GO:0009636', ' GO:0005739', ' GO:0004197', ' GO:0006508', ' GO:0003690', ' GO:0003697', 'GO:0000209', ' GO:0003729', ' GO:0042802', ' GO:0004177', ' GO:0000978', 'GO:0000122', ' GO:0046677', ' GO:0043418', ' GO:0008234', ' GO:0042493'}
>>> len(union)
21

# |intersection| / |union|
>>> float(len(intersection))/len(union)
0.23809523809523808
```

### PART B

The score of a protein s(p1, p2) where p1==p2 would be 1. Proof:

Within set theory, the idempotency law states that any union or intersection of any set with itself is equal to itself.
$A \cap A = A$
$A \cup A = A$

therefore:
$$s(A, A)=\frac{\left |A \cap A  \right |}{\left |A \cup A  \right |} = \frac{\left | A \right |}{\left | A \right |} = 1$$

## Question 7

**Q:**
The file 'SCOP_selection.txt' contains around 200 Uniprot IDs. For each pair of proteins, we want to determine whether they are functionally similar using the above scoring function. 

When you exclude pairs consisting of the same protein, how many unique pairs are there?  

**Answer:**
To determine the number of unique combinations between two equal sets of proteins, the following formula applies:

(n*(n-1))/2

The total is divided by two because combinations would otherwise be included twice assuming that any combination (as opposed to a permutation) is seen as equal (eg: (a,b)=(b,a)).

In the case of n being 209, the total number of unique combinations with same-protein pairs being excluded is:

(209*208)/2 = 21736

So a total of 21736 pairs would be scored for.

# Q

1. **Combinaties** of permutaties? Voor vraag: combinaties, voor code: permutaties.
2. bereik similarity **inclusive** of exclusive range