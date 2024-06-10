These files contain code for producing simplicial set structures for spaces homotopy equivalent to n-dimensional complex projective space, using the algorithm described in "Triangulations of complex projective spaces" by Sergeraert. The f-vectors for these models for CP^n:

n=2:  1  0  2   3   3
n=3:  1  0  3  10  25   30    15
n=4:  1  0  4  22  97  255   390   315    105
n=5:  1  0  5  40 271 1197  3381  5975   6405   3780    945
n=6:  1  0  6  65 627 4162 18496 54789 107933 139230 112770 51975 10395

Kenzo: 
- https://www-fourier.ujf-grenoble.fr/~sergerar/Kenzo/
- https://github.com/gheber/kenzo

The results for CP^2, CP^3, and CP^4 have been saved in the corresponding text files. The file S4.txt includes the 4-sphere, just for testing purposes. These files can be processed by the function "simplicial_data_from_kenzo_output" in sage/topology/simplicial_set.py. To get a simplicial set structure for CP^n using Kenzo in sbcl, do the following. 

;;
;; Start Kenzo.
;;
(require :asdf)
(require :kenzo)
(in-package "CAT")
;;
;; Define K(Z,2).
;;
(setf kz2 (k-z 2))
;;
;; Define effective homology version of K(Z,2).
;;
(setf efhm-kz2 (efhm kz2))
;;
;; The previous command produces output of the form
;; [K153 Homotopy-Equivalence K13 <= K143 => K139]
;;
;; In the following, replace "139" with the right-hand number.
;; Replace "4" with the desired dimension: 2n if you're constructing
;; CP^n.
;;
;; That is, the point is to find the smallest subcomplex of K(Z,2)
;; which contains the given homology class. If you replace "4" with
;; "2n", this should give a complex homotopy equivalent to CP^n.
;;
(chcm-homology-gen (k 139) 4)
(setf g (first *))
(setf z4 (lf efhm-kz2 (rg efhm-kz2 g)))
(multiple-value-setq (ssz4 incl) (gmsms-subsmst kz2 z4))
;;
;; Now ssz4 is a model for the 4-dimensional complex CP^2, so display
;; its nondegenerate simplices through dimension 4.
;;
(show-structure ssz4 4)
