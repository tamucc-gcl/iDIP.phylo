# iDIP.phylo
Function from Gaggiotti et al. (2018) Supplement 1 that works with HierDpart package

## Install Dependencies

```r
install.packages("ade4")
install.packages("ape")
install.packages("phytools")
install.packages("promises") 
install.packages("devtools")
library(devtools)
install_github("xinghuq/HierDpart")

library(ade4)
library(ape)
library(phytools)
library(HierDpart)
source(INSERT URL TO R SCRIPT)
```

## Creating matrix for `struct` argument

You need a matrix that follows the instructions in supplement 1, here's an example with a lot of levels

```r
> # Level 1: Regions (first 10 samples are Region1, next 10 are Region2)
> Region <- c(rep("Region1", 10), rep("Region2", 10))
> # Level 2: Archipelagos 
> # Region1: Arch1 for first 5 samples, Arch2 for next 5; Region2: Arch3 for first 5 samples, Arch4 for next 5.
> Archipelago <- c(rep("Arch1", 5), rep("Arch2", 5), rep("Arch3", 5), rep("Arch4", 5))
> # Level 3: Island Groups 
> # In each Archipelago, assign 2 groups: the first group gets 3 samples, the second gets 2.
> IslandGroup <- c(rep("IG1_1", 3), rep("IG1_2", 2),    # Region1, Arch1
+                  rep("IG2_1", 3), rep("IG2_2", 2),    # Region1, Arch2
+                  rep("IG3_1", 3), rep("IG3_2", 2),    # Region2, Arch3
+                  rep("IG4_1", 3), rep("IG4_2", 2))    # Region2, Arch4
> # Level 4: Islands
> # Each Island Group corresponds to exactly one Island (names match the group, but with an "I_" prefix).
> Island <- c(rep("I_1_1_1", 3), rep("I_1_1_2", 2),   # For IG1_1 and IG1_2 (Region1, Arch1)
+             rep("I_1_2_1", 3), rep("I_1_2_2", 2),   # For IG2_1 and IG2_2 (Region1, Arch2)
+             rep("I_2_1_1", 3), rep("I_2_1_2", 2),   # For IG3_1 and IG3_2 (Region2, Arch3)
+             rep("I_2_2_1", 3), rep("I_2_2_2", 2))   # For IG4_1 and IG4_2 (Region2, Arch4)
> # Level 5: Locations
> # To avoid having the same location name in two different islands, we incorporate the island name.
> # For islands with 3 samples, use _Loc1, _Loc2, _Loc3; for islands with 2 samples, use _Loc1 and _Loc2.
> Location <- c(
+   paste0("I_1_1_1_Loc", 1:3),    # 3 samples for I_1_1_1
+   paste0("I_1_1_2_Loc", 1:2),      # 2 samples for I_1_1_2
+   paste0("I_1_2_1_Loc", 1:3),      # 3 samples for I_1_2_1
+   paste0("I_1_2_2_Loc", 1:2),      # 2 samples for I_1_2_2
+   paste0("I_2_1_1_Loc", 1:3),      # 3 samples for I_2_1_1
+   paste0("I_2_1_2_Loc", 1:2),      # 2 samples for I_2_1_2
+   paste0("I_2_2_1_Loc", 1:3),      # 3 samples for I_2_2_1
+   paste0("I_2_2_2_Loc", 1:2)       # 2 samples for I_2_2_2
+ )
> # Level 6: Sample identifiers (each unique)
> Sample <- paste0("S", 1:20)
> # Combine all levels into a matrix where each column represents one sample
> Struc <- rbind(Region, Archipelago, IslandGroup, Island, Location, Sample)
> # View the structure matrix
> print(Struc)
            [,1]           [,2]           [,3]           [,4]          
Region      "Region1"      "Region1"      "Region1"      "Region1"     
Archipelago "Arch1"        "Arch1"        "Arch1"        "Arch1"       
IslandGroup "IG1_1"        "IG1_1"        "IG1_1"        "IG1_2"       
Island      "I_1_1_1"      "I_1_1_1"      "I_1_1_1"      "I_1_1_2"     
Location    "I_1_1_1_Loc1" "I_1_1_1_Loc2" "I_1_1_1_Loc3" "I_1_1_2_Loc1"
Sample      "S1"           "S2"           "S3"           "S4"          
            [,5]           [,6]           [,7]           [,8]          
Region      "Region1"      "Region1"      "Region1"      "Region1"     
Archipelago "Arch1"        "Arch2"        "Arch2"        "Arch2"       
IslandGroup "IG1_2"        "IG2_1"        "IG2_1"        "IG2_1"       
Island      "I_1_1_2"      "I_1_2_1"      "I_1_2_1"      "I_1_2_1"     
Location    "I_1_1_2_Loc2" "I_1_2_1_Loc1" "I_1_2_1_Loc2" "I_1_2_1_Loc3"
Sample      "S5"           "S6"           "S7"           "S8"          
            [,9]           [,10]          [,11]          [,12]         
Region      "Region1"      "Region1"      "Region2"      "Region2"     
Archipelago "Arch2"        "Arch2"        "Arch3"        "Arch3"       
IslandGroup "IG2_2"        "IG2_2"        "IG3_1"        "IG3_1"       
Island      "I_1_2_2"      "I_1_2_2"      "I_2_1_1"      "I_2_1_1"     
Location    "I_1_2_2_Loc1" "I_1_2_2_Loc2" "I_2_1_1_Loc1" "I_2_1_1_Loc2"
Sample      "S9"           "S10"          "S11"          "S12"         
            [,13]          [,14]          [,15]          [,16]         
Region      "Region2"      "Region2"      "Region2"      "Region2"     
Archipelago "Arch3"        "Arch3"        "Arch3"        "Arch4"       
IslandGroup "IG3_1"        "IG3_2"        "IG3_2"        "IG4_1"       
Island      "I_2_1_1"      "I_2_1_2"      "I_2_1_2"      "I_2_2_1"     
Location    "I_2_1_1_Loc3" "I_2_1_2_Loc1" "I_2_1_2_Loc2" "I_2_2_1_Loc1"
Sample      "S13"          "S14"          "S15"          "S16"         
            [,17]          [,18]          [,19]          [,20]         
Region      "Region2"      "Region2"      "Region2"      "Region2"     
Archipelago "Arch4"        "Arch4"        "Arch4"        "Arch4"       
IslandGroup "IG4_1"        "IG4_1"        "IG4_2"        "IG4_2"       
Island      "I_2_2_1"      "I_2_2_1"      "I_2_2_2"      "I_2_2_2"     
Location    "I_2_2_1_Loc2" "I_2_2_1_Loc3" "I_2_2_2_Loc1" "I_2_2_2_Loc2"
Sample      "S17"          "S18"          "S19"          "S20"
```

## Sources

* [Gaggiotti et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/eva.12593)
  * Supplements are the `*.pdf` in this repo
* [HierDpart](https://github.com/xinghuq/HierDpart)
