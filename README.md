# Goal

Goal is to benchmark KSW2 vs WFA2 on realistic PacBio CLR data. We mimic one
alignment step in the CCS algorithm, here the [draft
step](https://ccs.how/how-does-ccs-work.html#2-draft-generation). For this, all
subreads from one ZMW are aligned against a representative subread called
backbone. We perform mapping / seeding and extract alignment regions. The
histogram of those alignment regions has a sharp peak at ~250 bp, as we require
at least 200 bp long regions, but some are much longer due to the nature of
single-molecule sequencing.
