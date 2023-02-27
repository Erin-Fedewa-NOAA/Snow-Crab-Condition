# Snow-Crab-Condition
Script, data and output for NPRB funded project: Bering Sea snow crab lipid condition metrics in relation to temperature and recruitment dynamics. The "condition_haul_master.csv" contains immature snow crab sampled for hepatopancreas on 2019, 2021 and 2022 EBS and NBS NOAA bottom trawl surveys.

# Objectives:
1)	Determine spatial and interannual variation in body condition (i.e. total lipid content) and energy allocation (i.e. lipid class composition) of juvenile snow crab in the Bering Sea 
2)	Relate annual variation in lipid condition metrics to Bering Sea temperatures and annual estimates of snow crab recruitment in order to understand the importance of an energetic-based mechanism in linking thermal conditions to recruitment
3)	Assess the performance of a hepatosomatic index in predicting lipid content of field-collected snow crab for the development of a practical snow crab energetic condition index 

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
| End   | 2022-08-31 |

 
 

### Spatial

| LME     |                     |
| ------------- | ------------- |
|                | Eastern Bering Sea |
|               |   Northern Bering Sea |

![Rplot](https://user-images.githubusercontent.com/59858752/216689865-6766501b-98e4-420e-98cd-40a75e7a4ca1.png)


Attributes
----------
Three master datasets have been produced for further analyses via the "append haul" script. 1) total_FA_master.csv includes all snow crab biometric, fatty acid and haul level data, with data attributes listed below 2) percWT_FA_master.csv includes individual fatty acids for each hepatopancreas sample by % weight. Haul and biometric data in the attribute table below are included as well 3) perWWT_FA_master.csv includes individual fatty acid data reported as fatty acid per wet weight (micrograms/mg)

| Name    |    Description   |   Unit    |
| ------- | ---------------- | ---------- |
| `cruise` | Cruise ID for Bering Sea bottom trawl surveys. YYYY-01 indicates EBS surveys, YYYY-02 indicates NBS surveys. See RACEBASE or AKFIN for additional NOAA cruise metadata  |   ID Code
| `gis_station`   | Alpha-numeric designation for the station established in the design of AFSC standardized surveys | numeric/text                                      | `area_swept`   |   Unit of effort for AFSC bottom trawl surveys: computed by distance towed*mean net width   | numeric, in ha
|  `cpue`   |   Station-level snow crab density, calculated as CPUE"   |   numeric, crab/nmi^2
|  `vial_id`    |    Unique ID for hepatopancreas sample. AKK=Alaska Knight samples, V=Vesteraalen samples   |   numeric/text
|  `vessel`  |     ID number of the vessel used to collect data for that haul associated with vessel name    |   numeric
|  `haul`      |  Uniquely identifies a sampling event (haul) within an AFSC cruise. It is a sequential number, in chronological order of occurrence |  numeric
|  `sex` | Sex of specimen sampled. 1=Male, 2-Female   |  numeric
|  `cw`  |  Carapace width of specimen sampled | numeric, in mm
|  `ch_cc`   |   Chela height (males) or clutch code (females) used to determine maturity. 0 clutch code for females signifies 000, or immature   | numeric, in mm
 | `crab_wgt`  |   Whole crab weight of specimen sampled (prior to dissection). Only individuals with no missing limbs weighed | numeric, in g
|  `hepatotray_wt`  | 2019 samples only: Total weight of aluminum weigh boat    | numeric, in g
|  `hepato_wwt`    | 2019 samples only: Total weight of aluminum weigh boat PLUS 1-2 grams of hepatopancreas  | numeric, in g
|  `hepato_dwt`    |2019 samples only: Total weight of aluminum weigh boat PLUS 1-2 grams of hepatopancreas after drying at 70C for 48hrs   | numeric, in g
|  `muscletray_wt`  | 2019 samples only: Total weight of aluminum weigh boat | numeric, in g
|  `muscle_wwt`  |  2019 samples only: Total weight of aluminum weigh boat PLUS 1-2 grams of muscle tissue from 1st walking leg during dissection |  numeric, in g
|  `muscle_dwt`    |   2019 samples only: Total weight of aluminum weigh boat PLUS 1-2 grams of muscle tissue after drying at 70C for 48hrs   | numeric, in g
|  `process_date`  |   2019 samples only: date of dissection at KFRC (crab frozen whole in 2019 and dissected post-freeze)  | date
|  `notes`  |    Additional comments recorded during at-sea dissections | text
|  `maturity`  |  Maturity of specimen sampled. 0=Immature, 1=Mature"  | numeric                                                                                        | 'mid_latitude'       |   Latitude of specimen collection. Designates latitude at start of haul for AFSC standardized surveys    | numeric
|  `mid_longitude`    | Longitude of specimen collection. Designates longitude at start of haul for AFSC standardized surveys | decimal degree
 | `bottom_depth`    |    Bottom depth at station for AFSC standardized surveys  | numeric, in m
|  `gear_temperature`   |    Bottom temperature at sampling station | degree C
 | `year`     |        Year of specimen collection | numeric
|  `region`        |   Pre-defined spatial sampling strata in the Bering Sea. Regions 1-6 correspond with EBS BSIERP regions, Regions 7-9 correspond with NBS BSIERP regions  | numeric
|  `lme`         |   Large marine ecosystem: EBS=eastern Bering Sea, NBS=northern Bering Sea   |  text

:::

Distribution
------------

  File                                              Format    
  ------------------------------------------------- -------- -------------------------------------------------
  `Snow Crab Condition`   `csv`    [Download](https://github.com/Erin-Fedewa-NOAA)
