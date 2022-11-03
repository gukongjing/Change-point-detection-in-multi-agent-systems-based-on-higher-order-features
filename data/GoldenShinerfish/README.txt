This documentation file was generated on 2020-12-20 by Mathew Titus.

-------------------
# GENERAL INFORMATION
-------------------

1. Title of Dataset 

Fish Schooling Data Subset

2. Creator Information

Name: Yael Katz

Name: Kolbjørn Tunstrøm
Email: kolbjorn@chalmers.se

Name: Christos C. Ioannou

Name: Cristián Huepe

Name: Iain D. Couzin
Email: icouzin@princeton.edu

3. Contact Information 

Name: Mathew Titus
Institution: The Prediction Lab, LLC
Email: mat@thepredictionlab.com

-------------------
CONTEXTUAL INFORMATION
-------------------

1. Abstract for the dataset 

The data is a JSON format file containing the position, velocity, and fish identifier data for 300 golden shiners in a shallow (depth of 4.5 to 5 cm) rectangular water tank (2.1 by 1.2 meters). There are 5000 individual frames (samples of position and velocity) corresponding to video taken at a rate of 30 frames/s and analysed to extract individual fish’s trajectories. The fields px, py, vx, vy correspond to the x- and y-components of each detected fish’s position (at center of mass) and velocity. The onfish field gives each fish in the frame a  unique identifier. When an individual can no longer be distinguished, its identifier is retired; once the fish is again being tracked it is assigned a new ID. This data is a subsample of the frames acquired by the creators in the study https://doi.org/10.1073/pnas.1107583108, refer to it and its supplementary information for more details.

--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 

1. Licenses/restrictions placed on the data:

This work is licensed under the Creative Commons Attribution 4.0 International (CC-BY 4.0) license.

2. Links to publications related to the dataset:

Yael Katz, Kolbjørn Tunstrøm, Christos C. Ioannou, Cristián Huepe, and Iain D. Couzin. “Inferring the structure and dynamics of interactions in schooling fish.” Proceedings of the National Academy of Sciences. (2011) 108(46): 18720-18725. DOI: 10.1073/pnas.1107583108

Kolbjørn Tunstrøm, Yael Katz, Christos C. Ioannou, Cristián Huepe, Matthew J. Lutz, and Iain D. Couzin. “Collective states, multi stability and transitional behavior in schooling fish.” PLoS Computational Biology. (2013) DOI: 10.1371/journal.pcbi.1002915

Mathew Titus, George Hagstrom, and James R. Watson. “Unsupervised manifold learning of collective behavior.” PLoS Computational Biology. In review.

3. Recommended citation for the data:

Katz, Y., Tunstrøm, K., Ioannou, C.C., Huepe, C., and Couzin, I.D. (2021) Fish Schooling Data Subset [Dataset]. Oregon State University. https://doi.org/10.7267/zk51vq07c

4. Dataset Digital Object Identifier (DOI)

https://doi.org/10.7267/zk51vq07c

--------------------------
VERSIONING AND PROVENANCE
-------------------------- 

1. Last modification date

2021-01-05

2. Was data derived from another source?

Yes. The data was converted from a subset of a Matlab struct datafile, provided by K. Tunstrøm.

---------------------
DATA & FILE OVERVIEW
---------------------

1. File List
   A. Filename: schooling_frames.json
      Short description: A file in JSON-format with integer attributes (1, 2, …, 5000), each corresponding to a captured frame in chronological order. The value associated to each integer has attributes px, py, vx, vy, and onfish. The values of these attributes are arrays describing the fishes’ positions, velocities, and IDs, as described in the abstract above. 
