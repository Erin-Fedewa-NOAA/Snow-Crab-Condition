# Snow-Crab-Condition
Data and analyses for NPRB funded project #1911 "Bering Sea snow crab lipid condition metrics in relation to temperature and recruitment dynamics". This dataset includes immature snow crab sampled for hepatopancreas on 2019, 2021, 2022 and 2023 EBS and NBS NOAA bottom trawl surveys. Total fatty acid concentrations estimated from hepatopancreas samples were used to explore spatiotemporal variation in energetic condition in relation to temperature and snow crab density. 

# Project Objectives:
 1) Evaluate the performance of an indirect metric (hepatopancreas water content) in accurately estimating energetic condition of juvenile snow crab
 2) Assess spatial and temporal variation in energetic condition mid- and post-collapse
 3) Investigate the relative importance of bottom temperature and snow crab density on energetic condition of juvenile snow crab

# METADATA:
================================================

Bibliographic
-------------

| Published     | 02/02/2023   |
| ------------- | ------------- |
| Keywords      | snow crab |
|               |   Bering Sea |
|               |  fatty acids |
|               | condition |

 

Coverage
--------

### Temporal

| Begin    | 2019-06-01 |
| ------------- | ------|
| End   | 2023-08-31 |

 
 

### Spatial

| LME     |                     |
| ------------- | ------------- |
|                | Eastern Bering Sea |
|               |   Northern Bering Sea |

![Rplot](https://github.com/Erin-Fedewa-NOAA/Snow-Crab-Condition/blob/main/figures/data%20exploration/n_year.png?raw=true)



Attributes
----------
One master datasets has been produced for further modeling via the "append haul" script. "total_FA_master.csv" includes all snow crab biometric, total fatty acid and haul level data, with data attributes listed below. 

| Name    |    Description   |   Unit    |
| ------- | ---------------- | ---------- |
| `year`     |        Year of specimen collection | numeric
|  `Total_FA_Conc_WWT` |   total fatty acid wet weight concentration of hepatopancreas sample   |  mg FA/ g WWT
|  `vial_id`    |    Unique ID for hepatopancreas sample by year. AKK=Alaska Knight samples, V=Vesteraalen samples   |   numeric/text
| `cruise` | Cruise ID for Bering Sea bottom trawl surveys. YYYY-01 indicates EBS surveys, YYYY-02 indicates NBS surveys. See RACEBASE or AKFIN for additional NOAA cruise metadata  |   ID Code
| `gis_station`   | Alpha-numeric designation for the station established in the design of AFSC standardized surveys | numeric/text   
|  `total_benthic_cpue`    |    Average density of benthic invertebrates on EBS bottom trawl survey, filtered to know prey items of snow crab   |   numeric, in thous metric tons
| `area_swept`   |   Unit of effort for AFSC bottom trawl surveys: computed by distance towed*mean net width   | numeric, in ha
|  `cpue`   |   Station-level snow crab density, calculated as CPUE"   |   numeric, crab/nmi^2
|  `vessel`  |     ID number of the vessel used to collect data for that haul associated with vessel name    |   numeric
|  `haul`      |  Uniquely identifies a sampling event (haul) within an AFSC cruise. It is a sequential number, in chronological order of occurrence |  numeric
|  `sex` | Sex of specimen sampled. 1=Male, 2-Female   |  numeric
|  `cw`  |  Carapace width of specimen sampled | numeric, in mm
|  `ch_cc`   |   Chela height (males) or clutch code (females) used to determine maturity. 0 clutch code for females signifies 000, or immature   | numeric, in mm
 | `crab_wgt`  |   Whole crab weight of specimen sampled (prior to dissection). Only individuals with no missing limbs weighed | numeric, in g
|  `hepato_wwt`    | Total wet weight of hepatopancreas (subtracted from weight of weigh boat)  | numeric, in g
|  `hepato_dwt`    |Total dry weight of hepatopancreas after drying at 70C for 48hrs (subtracted from weight of weigh boat)  | numeric, in g
|  `maturity`  |  Maturity of specimen sampled. 0=Immature, 1=Mature"  | numeric                                                                                        
| 'mid_latitude'       |   Latitude of specimen collection. Designates latitude at start of haul for AFSC standardized surveys    | numeric
|  `mid_longitude`    | Longitude of specimen collection. Designates longitude at start of haul for AFSC standardized surveys | decimal degree
 | `bottom_depth`    |    Bottom depth at station for AFSC standardized surveys  | numeric, in m
|  `gear_temperature`   |    Bottom temperature at sampling station | degree C
 | `start_date`     |        Date of sampling | date, month/day/year
|  `DWT_WWT` |   hepatopancreas dry weight/hepatopancreas wet weight   |  numeric, in g
|  `Perc_DWT` |   hepatopancreas dry weight/ hepatopancreas wet weight x 100; percentage hepatopancreas dry weight    |  percentage
|  `Total_FA_Conc_DWT` |   total fatty acid dry weight concentration of hepatopancreas sample   |  mg FA/ g DWT
|  `WWT_DWT` |   hepatopancreas wet weight/dry weight    |  numeric, in g
|  `lme` |   Large marine ecosystem: EBS=eastern Bering Sea, NBS=northern Bering Sea   |  text

:::

Distribution
------------

  File                                              Format    
  ------------------------------------------------- -------- -------------------------------------------------
  `total_FA_master`   `csv`    [Download](https://github.com/Erin-Fedewa-NOAA)
