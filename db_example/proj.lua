
-- ++++++++++++++++++++++++++++++++++++++++ Model Data +++++++++++++++++++++++++++++++++++++++++++++++++++

model = getdata("Input JBOOST WAK-DP-A-LILA-MySE18-260-v0.lua")
model()
wave = (getdata("wave.lua"))
wave()
wind = getdata("wind.lua")
wind()

-- ++++++++++++++++++++++++++++++++++++++++ Configurations +++++++++++++++++++++++++++++++++++++++++++++++++++

os_LoadConfigGeneral(
{
	Config_name	= "C0_I",
	Result_name	= "Result DP-A",
	Scatter_name	= "EMPTY",
	FeModel_name	= "WAK",
	Wind_name	= "WindLoc1",
	Ocean_name	= "test",
}
)

-- ++++++++++++++++++++++++++++++++++++++++ FEModul Konfiguration ++++++++++++++++++++++++++++++++++++++++++++++++++++++
os_LoadConfigFEModul(
{
	foundation_superelement	= 501,	-- Knotennummer (-1=keine Ersatzmatizen fuer Pfahl, >0 Ersatzsteifigkeit bzw. -massen)

	found_stiff_trans	=  4.84E+09,
	found_stiff_rotat	=  7.46E+11,
	found_stiff_coupl	=  -4.16E+10,
	found_mass_trans	=  0,
	found_mass_rotat	=  0,
	found_mass_coupl	=  0,

	shearCorr	= 2.0,	-- Schubkorrekturfaktor (2.0 fuer Kreisquerschnitt und Timoshenko 0.0 fuer Bernoulli)
	res_NumEF	= 10,	-- Anzahl der Eigenfrequenzen, die ausgegeben werden
	hydro_add_mass	= 1,	-- Added Mass fuer analysen der Strukturdynamik 1 = mit added mass
}
)
os_LoadConfigOcean(
{
	water_density 	= 1027,	-- water density in [kg/m ]
	water_level 	= 0.86,	-- wrt LAT in [m] !!! there must be a node @ this height !!!
	seabed_level 	= -40.3,	-- wrt LAT in [m]
	growth_density 	= 1325,	-- marine growth density in [kg/m ]
}
)

os_OceanMG{id=HD, topMGsection = 2.86,   bottomMGsection = -9.14,    t=0.150} -- marine growth layer spec.
os_OceanMG{id=HD, topMGsection = -9.14,   bottomMGsection = -40.3,    t=0.100} -- marine growth layer spec.

os_LoadConfigFrequModul(
{
	frequRange		= 2,	-- 0 Hz to frequRange in Hz
	frequRangeSolution	= 300,	-- no. of frequency components
	maxEFcalc		= 10,	-- max number of considered frequencies
	damping_struct		= 0.0075,	-- damping ratio structural (material, viscous, soil)
	damping_tower		= 0.01,
	damping_aerodyn 	= 0.12,	-- damping ratio aerodynamic (weighted mean value over wind speeds)
	design_life		= 31.0,	-- structural design life in years
	N_ref			= 1.0E7,	-- no. of cycles for DEL calculation
	tech_availability	= 0.92,	-- technical availability turbine in percent --> 0.95*30y/31.0y
	-- Note that this coers also the commissioning and decommissioning times
	SN_slope		= 4,	-- S-N slope material
	h_refwindspeed		= 156.6,	-- height reference in m of wind speed data in wind wave scatter
	h_hub			= 156.6,	-- hub height wrt LAT in m
	height_exp		= 0.11,	-- height exponent
	-- height exp was set specifically here so that the average wind speed from the ROUGH Vw Hs scatter is 10.5 m/s at hub height
	TM02_period		= 0,	-- [0] if peak periods stated, [1] if zero crossing Tm02 periods stated -> noch nicht eingebaut...FOs
	refineScatter		= 0,	-- refinement of scatter data on Tp axis (factor 10 coded)  [0]=off; [1]=on; [2]=on incl. validation plots!
	fullScatter		= 1, 	-- calc full real scatter FLS for validation purpose
	WindspecDataBase 	= [[-]],	-- path on windspec database (*.csv) or "-" in case no wind load is to be considered
	res_Nodes		= {501.0},
}
)

-- ----------------------------------------------
-- Exexute program steps
-- ----------------------------------------------

--os_RunFEModul()
--os_RunFrequencyDomain()
os_RunHindcastValidation()

-- Ausgabe der Textdateien
--os_WriteResultsText([[Results_OLFD_Text]])
--os_PlotResultsGraph([[Result_OLFD_Graph]])
os_WriteResultsText([[Results_JBOOST_Text]])
os_PlotResultsGraph([[Result_JBOOST_Graph]])

