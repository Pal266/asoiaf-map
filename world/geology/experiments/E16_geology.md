# E16 — Rock ages, faults and seismic hazard

**Goal:** three things the resource layer needs, built from the plates (E06), the geological history (E07) and the
lithology (E14):
- the formation age of the surface rock of every land cell;
- the faults, with their type and the age of their last movement;
- a seismic hazard map for the game.

## Research

- **Ages of exposed rock.** On global geological maps the area of exposed rock per million years falls with age as
  a power law. The exponent is −0.86 on the FAO/UNESCO map and −0.85 on the GSC map (Wilkinson et al. 2009,
  "Global geologic maps are tectonic speedometers"). Integrated from 0 to 4000 Ma, this gives these shares of the
  land: Cenozoic 56 %, Mesozoic 12 %, Paleozoic 8 %, Precambrian 25 %. For comparison, the Cenozoic covers nearly
  half of the conterminous United States (USGS, Professional Paper 904).
- **Periods:** International Chronostratigraphic Chart (ICS, v2023/09), with the Precambrian by eon and era.
- **Plate boundaries and earthquakes.** PB2002 (Bird 2003) classifies the boundaries into seven step classes:
  subduction, oceanic and continental convergent, oceanic and continental transform, continental rift, oceanic
  spreading ridge. It also maps 13 orogens, the zones of distributed active deformation. The length-weighted
  median relative velocities are: subduction 64, spreading 34, oceanic transform 32, continental convergent 21,
  continental transform 16, continental rift 11, oceanic convergent 9 mm/yr. The seismic moment released on a
  fault grows with its slip rate.
- **Seismic hazard reference:** the GSHAP map (Giardini et al. 1999), peak ground acceleration with 10 %
  probability of exceedance in 50 years, in four classes: low < 0.08 g, moderate 0.08–0.24, high 0.24–0.40, very
  high > 0.40 g.

## Method (scripts/geology.py)

### Rock ages

Each surface rock of the lithology layer is dated by its setting. The setting comes from the province and event age
the lithology drew it from: the history read through the lithology's warp, which the lithology now saves.
- **Crystalline rocks** (metamorphic, plutonic, shield volcanics) take the age of their province's last tectonic
  event, up to 10 % younger (later intrusions).
- **Sediments** are sampled from the power law (exponent 0.86) between the limits of their setting:
  - platform cover between its basement's age and 2.58 Ma;
  - in old and young orogens, folded sediments up to 200 Myr older than the orogeny; molasse and cover younger than
    it (half of the sediments in old orogens, 30 % in young ones);
  - arc basins 2.58–100 Ma; plateau sediments older than the uplift; basin fill 0–66 Ma.
- **Volcanics:**
  - arc and rift volcanics are as young as their arc or rift;
  - flood basalts on platforms follow the cover law;
  - the Doom's pyroclastics in the Fourteen Flames are ~400 years old (canon).
- **Unconsolidated sediments and surface deposits:** unconsolidated sediments are Pliocene to Holocene, deposits
  Quaternary.

One smooth random field supplies the quantiles, so ages form coherent outcrop belts.

### Faults

Faults are drawn as one-cell lines, each with a type, the age of its last movement and an activity class:
- **Plate contacts** from the plates layer, typed by the boundary design:
  - trench: subduction megathrust;
  - ridge: spreading ridge in oceanic crust, normal fault (continental rift) in continental crust;
  - transform.

  All are active.
- **Contacts of provinces on continental crust**, on the younger side:
  - young orogens, arcs and plateaus: thrusts;
  - rifts and basins: normal faults;
  - old orogens against cratons and platforms: old thrust fronts;
  - other contacts between blocks of different age: sutures.

  The fault takes the age of the younger side's event. It is active if that is within 30 Myr, dormant within
  66 Myr, and inactive if older. Faults of rifted margins (the Narrow Sea, 25 Ma) went dormant when spreading
  began. The post-glacial canon features (the Neck, the Stepstones) are not tectonic and get no faults.

### Seismic hazard (scripts/hazard.py)

The index is the strongest of these contributions:
- each plate boundary: w_class × min(1, v / 20 mm/yr)^0.5 × exp(−d / 500 km). The weights are 1.0 for
  subduction, convergent and continental transform boundaries, 0.8 for continental rifts and 0.5 for oceanic
  ridges and transforms;
- active deforming zones, weight 0.5 inside and falling off outside. They are the young orogens, arcs and plateaus
  (≤ 30 Ma) and the continental rifts (1–15 Ma);
- dormant zones (30–66 Ma), 0.35 inside.

Velocities on world boundaries:
- ridges use their spreading rate from the history;
- trenches and transforms, whose rate the scenario does not give, take the PB2002 median of their class.

The GSHAP classes are cut at 0.30, 0.45 and 0.55. Rifted passive margins are not deforming zones: on Earth, PB2002
counts no such margin as an orogen.

## Earth test of the hazard (scripts/earth_hazard_validation.py)

The test uses the PB2002 steps and orogens and 71 cities with their GSHAP class. The classes were read from the
published map, and two classes are accepted near a class boundary.

| Run | Change | Cities right |
|---|---|---|
| a | distance to boundaries by class, decay 100 km | 0.58 |
| b | orogen polygons across the antimeridian unwrapped (New Hebrides–Fiji had covered all longitudes, Brasília included); velocity weighting (the Somalia–Africa boundary at ~1 mm/yr had made Johannesburg very high) | 0.65 |
| c | class thresholds searched, all weights free | 0.79 (leave-one-out 0.63: overfitted) |
| **d** | **boundary weights fixed by class; search only the decay length, reference velocity, zone weight and thresholds** | **0.87 (leave-one-out 0.85)** |

Misses in run d:
- underestimated: Anchorage, Almaty and Bishkek (very high, modelled high); Nairobi and Vladivostok (moderate,
  modelled low);
- overestimated: Athens, Delhi, Karachi and Shanghai (one class too high).

The PB2002 data come from the public copy at github.com/fraxen/tectonicplates (ODC-BY). They are read for the test
only and are not stored in the repository.

## Rock ages on the world (run d)

| Era | World | Power law (Wilkinson et al. 2009) |
|---|---|---|
| Cenozoic | 0.559 | 0.563 |
| Mesozoic | 0.093 | 0.116 |
| Paleozoic | 0.119 | 0.076 |
| Precambrian | 0.229 | 0.245 |

- **Area-age slope**, from area per Myr in log bins from 1 to 3000 Ma: −0.97 (Earth −0.86). The world is a little
  steeper because the crystalline ages cluster at the scenario's events (300, 420, 600 and 1800 Ma). The 600 Ma
  mobile belts of the cratons alone make the Neoproterozoic 16 % of the land.
- **The Cenozoic and Precambrian shares are not a free test.** Sediments are drawn from the same law. What the
  comparison does check is the mix of settings: how much land is old basement, folded belts, cover and young
  deposits.
- **Periods on land:** Quaternary 27.6 %, Neogene 19.2 %, Paleogene 9.1 %, Cretaceous 5.5 %, Jurassic 2.3 %,
  Triassic 1.5 %, Permian 1.4 %, Carboniferous 5.8 %, Devonian 1.7 %, Silurian 0.6 %, Ordovician 1.1 %, Cambrian
  1.3 %, Neoproterozoic 16.4 %, Mesoproterozoic 1.5 %, Paleoproterozoic 4.7 %, Archean 0.3 %.
- **Westeros south of the Wall:** Quaternary and Neogene 52 % (deposits and young cover); Paleogene and Cretaceous
  19 %; Devonian, Carboniferous and Permian 13 % (the worn Caledonian and Variscan hills).
- **Beyond the Wall:** Neogene cover over Paleoproterozoic and Neoproterozoic shield.
- **The places:** Castle Black stands on Mesoproterozoic limestone (shield cover); Winterfell and King's Landing on
  Quaternary sediments; Sunspear on Neogene limestone; Asshai on 5 Ma basalt.

## Faults on the world (run d)

| Type | Active, km | Dormant, km | Inactive, km |
|---|---|---|---|
| subduction megathrust | 49 800 | | |
| spreading ridge | 51 800 | | |
| transform | 6 500 | | |
| thrust | 98 100 | 2 600 | 177 200 |
| normal | 16 500 | 30 600 | |
| suture | | | 15 700 |

- **Westeros:** almost all its faults are inactive Caledonian and Variscan thrust fronts and sutures. The Red
  Mountains (40 Ma) are dormant, and the faults of its east coast went dormant with the opening of the Narrow Sea.
- **Essos:** active thrusts in the southern collision belt (Painted Mountains, Ghiscar, Lhazar) and the Bone
  Mountains; normal faults around the Bleeding Sea rift.

## Seismic hazard on the world (run d)

| Region | Low | Moderate | High | Very high |
|---|---|---|---|---|
| All land | 70 % | 9 % | 19 % | 2 % |
| Westeros south of the Wall | 72 % | 18 % | 4 % | 5 % |
| Beyond the Wall | 93 % | 5 % | 1 % | 1 % |
| Western Essos | 42 % | 22 % | 33 % | 4 % |

- **Very high:**
  - southern Dorne and Sunspear, on the Dorne transform between Westeros and Sothoryos (like Gibraltar–Algeria);
  - Volantis and the Valyrian peninsula, above the Summer Sea trench;
  - Asshai.
- **High:** the southern collision belt of Essos.
- **Moderate:** the Narrow Sea coasts and the Red Mountains.
- **Low:** the North, the Riverlands, the Reach and King's Landing.

Run d is deterministic (two runs are identical) and takes ~110 s.

## Still open

- **Uniform ridge rates.** The plate design gives one rate per ridge and no poles of rotation, so the Narrow Sea
  ridge spreads at 20 mm/yr along its whole length. Toward its pole it should slow, which would lower the hazard
  near Braavos and in the far north.
- **Straight terrane edges.** The canon hills and mountains of the relief (E04) are hand-drawn polygons. As
  terranes of their own age (420 and 300 Ma) they keep fairly straight edges in the age map, even through the
  lithology's 0.8° warp.
- **No ore bodies yet.** They are the next layer, built on these ages, faults and the lithology.

## Addendum (2026-11-22, from E17): shear zones of the Westerlands

The Westerlands hills, an accreted arc terrane since decision 13, carry a dense network of shear zones: the zero
lines of a smooth random field, ~30 km apart, as inactive thrusts of the terrane's age (300 Ma). They add 6 300 km
of inactive thrusts (177 200 → 183 500 km). Hazard and ages are unchanged.

## Sources

- Wilkinson, B. H., McElroy, B. J., Kesler, S. E., Peters, S. E., Rothman, E. D. (2009), Global geologic maps are
  tectonic speedometers — Rates of rock cycling from area-age frequencies, GSA Bulletin 121, 760.
- USGS Professional Paper 904, The Cenozoic rocks: a discussion to accompany the Geologic Map of the United States.
- Cohen, K. M., Finney, S. C., Gibbard, P. L., Fan, J.-X., International Chronostratigraphic Chart v2023/09.
- Bird, P. (2003), An updated digital model of plate boundaries, G3 4(3), 1027; public copy
  github.com/fraxen/tectonicplates (H. Ahlenius, ODC-BY).
- Giardini, D., Grünthal, G., Shedlock, K. M., Zhang, P. (1999), The GSHAP Global Seismic Hazard Map, Annali di
  Geofisica 42, 1225.
