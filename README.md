# ribosome-building-analysis
Code used in the study: "Insights into Translocation Mechanism and Ribosome Evolution from Cryo-EM Structures of Translocation-Intermediates of Giardia intestinalis."

## svd_axis_angle
A set of scripts to determine the rotation characeristics of translocation intermediates in a pairwise fashion.
Measure the 3-dimensional distance between rRNA residues of two structures.
Determine the rotational axis and rotation vectors (between phosphates for instance) between two intermediate states.
Compare rotational axes between transitions.

![Pairwise Euclidean norm](./images/SI_eucnorm.png)
3-dimensional distance between rRNA residues from one translocation intermediate to another (left) with the small subunit (SSU) body residues aligned (right).

![Rotation axes comparison](./images/rotation_axes_01.png)
Comparison of ratcheting and rolling rotational axes for Giardia translocation intermediates.

## coot_progressive_refinement
A script to perform real space refinement in Coot using a docked structure. This script orders the residues to be refined by their distance from the (macro)molecule's centre of mass. Then it performs this refinement in a step-wise fashion in a sphere of a given radius (determined by the user, some choices in the Coot menu are 8 and 25 Å).
