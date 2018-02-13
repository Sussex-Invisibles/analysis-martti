/////////////////////////////////////////
// Original by Gwenaelle Lefeuvre
// g.lefeuvre@sussex.ac.uk
// Modified by Ingrida Semenec
// isemenec@laurentian.ca
// Modified by Martti Nirkko
// m.nirkko@sussex.ac.uk
/////////////////////////////////////////
//
// Geometry of AmBe Source "High"
// Source: J. Loach PhD thesis, 2008, p72
// 
// Position of the source: 0, 0, 17.5
//
// Delrin: material of the SNO source. NOT
// COMPATIBLE with LS!
//
// Envelope volume has no physical meaning,
// it simply serves a mother volume for the
// other components of the source. By default
// it is invisible.
//

/////////////////////////////////////////
// Part 0, "Envelope" volume
//
{
type: "GEO",
version: 1,
index: "AmBeSourceShielding",
run_range:  [0, 0],
pass: 0,
comment: "",
timestamp: "",
enable: 1,

mother: "inner_av",

factory: "solid",
solid: "polycone",

z:       [ 0.4, 60.7, 60.7, 73.8, 73.8, 173.5, 369.3, 369.3, 378.5],
r_inner: [ 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0],
r_outer: [22.2, 22.2, 14.8, 14.8,  7.2,  14.4,  14.4,  31.9,  31.9],

position: [0.0, 0.0, 0.0],

//material: "labppo_scintillator",    //same material as mother volume
material: "air",                    //same for space between source and lead
vis_color: [0.6, 0.3, 1.0, 0.2],
vis_invisible: 1,
}

/////////////////////////////////////////
// Part 1, Stainless Steel Outer Shielding
//
{
type:  "GEO",
version: 1,
index: "AmBeCan",
run_range: [0, 0],
pass: 0,
comment: "",
timestamp: "",
enable: 1,

mother: "AmBeSourceShielding",

factory: "solid",
solid: "polycone",

z:       [ 0.5,  2.5,  2.5, 57.6, 57.6, 60.6, 60.6, 66.95], //1mm + 2 mm for screw + 6.35 mm on top for connector
r_inner: [ 0.0,  0.0, 21.0, 21.0,  0.0,  0.0,  0.0,  0.0],
r_outer: [22.13, 22.13, 22.13, 22.13, 22.13, 22.13, 14.72, 14.72],

position: [0.0, 0.0, 0.0],        // wrt AmBeSourceShielding

material: "stainless_steel",      // delrin not compatible with LAB
vis_color: [0.4, 0.4, 0.4, 0.1],  // light gray
}

////////////////////////////////////////
// Part 2, Lining. Lead
{
type:  "GEO",
version: 1,
index: "AmBeShield",
run_range: [0, 0],
pass: 0,
comment: "",
timestamp: "",
enable: 1,

mother: "AmBeSourceShielding",

factory: "solid",
solid: "polycone",

z:       [ 2.5,  4.5,  4.5, 55.6, 55.6, 57.6], //2 mm
r_inner: [ 0.0,  0.0, 19.0, 19.0,  0.0,  0.0],
r_outer: [21.0, 21.0, 21.0, 21.0, 21.0, 21.0],

position: [0.0, 0.0, 0.0],        // wrt AmBeSourceShielding

material: "Lead",                 // lead or hevimet
vis_color: [0.6, 1.0, 1.0, 0.3],  // light teal
}

////////////////////////////////////////
// Part 3, Actual Source Body
{
type:  "GEO",
version: 1,
index: "AmBeBody",
run_range: [0,0],
pass: 0,
comment: "",
timestamp: "",
enable: 1,

mother: "AmBeSourceShielding",

factory: "solid",
solid: "polycone",

z:       [ 5.0, 39.8, 39.8, 49.8, 49.8, 55.1],
r_inner: [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
r_outer: [14.6, 14.6, 16.2, 16.2, 14.6, 14.6],

position: [0.0, 0.0, 0.0],        // wrt AmBeSourceShielding

material: "G4_POLYOXYMETHYLENE",  // delrin
vis_color: [0.8, 0.0, 1.0],       // pink

}

/////////////////////////////////////////
// Part 4, Teflon Stem
{
type:  "GEO",
version: 1,
index: "AmBeHolder",
run_range: [0, 0],
pass: 0,
comment: "",
timestamp: "",
enable: 1,

mother: "AmBeSourceShielding",

factory: "solid",
solid: "polycone",

z:       [67.35, 73.75, 73.75, 73.75, 173.45, 369.35, 369.35, 378.45],
r_inner: [  0.0,   0.0,   0.0,  5.05,   5.05,   5.05,    0.0,    0.0],
r_outer: [14.72, 14.72,  7.15,  7.15,  14.30,  14.30,  31.81,  31.81],

position: [0.0, 0.0, 0.0],        // wrt AmBeSourceShielding

material: "ptfe",                 // teflon
vis_color: [0.4, 0.4, 0.4, 0.1],  // light grey
}

/////////////////////////////////////////
// Part 5, Neutron absorber rod (optional)
{
type:  "GEO",
version: 1,
index: "AmBeAbsorber",
run_range: [0, 0],
pass: 0,
comment: "",
timestamp: "",
enable: 1,

mother: "AmBeSourceShielding",

factory: "solid",
solid: "polycone",

z:       [74.15, 368.95],
r_inner: [ 0.0,    0.0],
r_outer: [ 5.0,    5.0],

position: [0.0, 0.0, 0.0],        // wrt AmBeSourceShielding

material: "TEMPLATE"              // neutron absorber for high-energy gammas
                                  // use same material as stem to ignore!
vis_color: [1.0, 0.6, 0.0, 0.2],  // light orange
}
//
//
/////////////////////////////////////////
