# November, 2020
# Kevin Dorma
# module for common PSV calculations
# Dependancies?
# GridInterpolations
# next, HEM

module PSVrelief

using DataFrames
using GridInterpolations

export PSVvaporFlux, PSVvaporRate, PSVvaporSize, PSVliquidSize, PSVliquidRate
export thermExpansionRate, poolFireReliefRate, PSVsteamFlux, PSVsteamRate, PSVsteamSize, liquidVaporizeReliefRate
export waterPsat, waterTsat, getKsh, PSVareaOrifice, PSVfindOrifice

# to do
# provide either PSV descriptor or area for rating.
# create function to return the next larger PSV descriptor
# implement proper steam equations
# implement HEM, Homogeneous Equilibrium Method for steam relief

function waterTsat(P_kPa)
    # equation for saturated steam, 1/T = A + B*ln P + C/ln P
    # I fit this myself, 100 to 20000 kPaa
    # pressure is in kPaa
    # we need to replace this with IAPWS IFC97 formulas

    lnP = log(P_kPa)
 
    # constants for fit of 1/T = A + B*ln P + C/ln P
    tA = 0.00379302
    tB = -0.000220828
    tC = -0.000425693

    invT = tA + tB*lnP + tC/lnP
    return((1.0/invT) - 273.15)
end


    
    
function waterPsat(T_c)
    # equation for saturated steam, ln P_Pa = A + B/Tk + C*Tk + D*ln(Tk)
    # I fit this myself, 100 to 20000 kPaa
    # pressure is in kPaa

    T_k = T_c + 273.15

    pA = 116.6408494
    pB = -8572.035364
    pC = 0.013736471
    pD = -14.73621925
    lnPsat = pA + pB/T_k + pC*T_k + pD*log(T_k)
    return(exp(lnPsat)/1000.0)

end

function getKsh(PkPa, State)
    # State is either a temperature in degC or a string
    # if this is a temperature, then look up the value in the table
    # if not, just return 1.0 because this is saturated
    # throw an error if the tempeature is below saturation
    # tables for superheat factors. I had a larger table but it threw weird errors. Smaller seems to be better.
    
    Ksh_tC = fill(0.0,(10))
    Ksh_pkPa = fill(0.0,(26))
    Ksh_table = fill(0.0, (26,10))    # pressure in rows, temperature in columns

# julia code for Ksh temperature values in degC
    Ksh_tC = [93.33333333	,
		148.8888889	,
		204.4444444	,
		260.0	,
		315.5555556	,
		371.1111111	,
		426.6666667	,
		482.2222222	,
		537.7777778	,
		565.5555556	];

    # julia code for Ksh pressure values in kPa   
    Ksh_pkPa = [137.8951817	,
		344.7379543	,
		689.4759087	,
		1034.213863	,
		1378.951817	,
		1723.689772	,
		2068.427726	,
		2413.16568	,
		2757.903635	,
		3102.641589	,
		3447.379543	,
		3792.117498	,
		4136.855452	,
		4826.331361	,
		5515.807269	,
		6205.283178	,
		6894.759087	,
		7584.234995	,
		8273.710904	,
		8963.186813	,
		9652.662721	,
		10342.13863	,
		12065.8284	,
		13789.51817	,
		17236.89772	,
		20684.27726	];


    # Julia code for Ksh, rows are pressure, columns are temperature    
    # subcooled values are denoted with KSH = 1.
    # Interpolation near the saturation temperature could give artificially low value for Ksh
    # a better method might be to replace the value 1 for the nearest subcooled temperature
    # with the value that gives 1 when interpolated to the saturation tempeature.
    # Clumsy, but it should work
    Ksh_table = [1	0.99455814	0.987	0.93	0.882	0.841	0.805	0.774	0.745	0.732	;
		1	0.997925224	0.987	0.93	0.882	0.841	0.805	0.774	0.745	0.732	;
		1	1	0.998	0.935	0.885	0.843	0.807	0.775	0.746	0.733	;
		1	1	0.984	0.94	0.888	0.846	0.808	0.776	0.747	0.733	;
		1	1	0.979	0.945	0.892	0.848	0.81	0.777	0.748	0.734	;
		1	1	1	0.951	0.895	0.85	0.812	0.778	0.749	0.735	;
		1	1	1	0.957	0.898	0.852	0.813	0.78	0.75	0.736	;
		1	1	1	0.963	0.902	0.854	0.815	0.781	0.75	0.736	;
		1	1	1	0.963	0.906	0.857	0.816	0.782	0.751	0.737	;
		1	1	1	0.961	0.909	0.859	0.818	0.783	0.752	0.738	;
		1	1	1	0.961	0.914	0.862	0.82	0.784	0.753	0.739	;
		1	1	1	0.962	0.918	0.864	0.822	0.785	0.754	0.74	;
		1	1	1	0.964	0.922	0.867	0.823	0.787	0.755	0.74	;
		1	1	1	1	0.931	0.872	0.827	0.789	0.757	0.742	;
		1	1	1	1	0.942	0.878	0.83	0.792	0.759	0.744	;
		1	1	1	1	0.953	0.883	0.834	0.794	0.76	0.745	;
		1	1	1	1	0.959	0.89	0.838	0.797	0.762	0.747	;
		1	1	1	1	0.962	0.896	0.842	0.8	0.764	0.749	;
		1	1	1	1	0.966	0.903	0.846	0.802	0.766	0.75	;
		1	1	1	1	0.973	0.91	0.85	0.805	0.768	0.752	;
		1	1	1	1	0.982	0.918	0.854	0.808	0.77	0.754	;
		1	1	1	1	0.993	0.926	0.859	0.811	0.772	0.755	;
		1	1	1	1	1	0.94	0.862	0.81	0.77	0.752	;
		1	1	1	1	1	0.952	0.861	0.805	0.762	0.744	;
		1	1	1	1	1	0.951	0.852	0.787	0.74	0.721	;
		1	1	1	1	1	1	0.831	0.753	0.704	0.684	];

        # we will use linear interpolation with P and T. 
        # Using ln P and 1/T might be more robust and needs to be investigated.


    linear_Ksh = reshape(Ksh_table,(10*26,1));
    Ksh_grid = GridInterpolations.RectangleGrid(Ksh_pkPa, Ksh_tC);  	# rectangular grid

    
    Ksh = 1.0                # default value
    if (typeof(State)<:Number) # if the State is some kind of number, do the following
        Ksh = GridInterpolations.interpolate(Ksh_grid,linear_Ksh,[PkPa,State])
        # check if we are subcooled
        pSat = waterPsat(State)
        if (pSat < PkPa)
            error("Temperature is in subcooled region")
        end
    end
    return (Ksh)
end

function PSVareaOrifice(letter)
# psvTable[1,:areaMM2]
# from the PSV designation letter, output the PSV area in MM2
    psvTable = DataFrame(Designation = String[], Flanges = String[], areaIN2 = Float64[], areaMM2 = Float64[])
    push!(psvTable, ("D","1.5D2",0.11,70.9676))
    push!(psvTable, ("E","1.5E2",0.20,126.45136))
    push!(psvTable, ("F","1.5F3",0.31,198.06412))
    push!(psvTable, ("G","1.5G3",0.50,324.51548))
    push!(psvTable, ("H","2H3",0.79,506.4506))
    push!(psvTable, ("J","3J4",1.29,830.32092))
    push!(psvTable, ("K","3K4",1.84,1185.80408))
    push!(psvTable, ("L","4L6",2.85,1840.64148))
    push!(psvTable, ("M","4M6",3.60,2322.576))
    push!(psvTable, ("N","4N6",4.34,2799.9944))
    push!(psvTable, ("P","4P6",6.38,4116.1208))
    push!(psvTable, ("Q","6Q8",11.05,7129.018))
    push!(psvTable, ("R","6R10",16.00,10322.56))
    push!(psvTable, ("T","8T10",26.00,16774.16))

    return(psvTable[psvTable[:,:Designation] .== letter, :areaMM2])
    
end

function PSVfindOrifice(area)
    # from the required area (mm2), find the index for the next size up
    psvTable = DataFrame(Designation = String[], Flanges = String[], areaIN2 = Float64[], areaMM2 = Float64[])
    push!(psvTable, ("D","1.5D2",0.11,70.9676))
    push!(psvTable, ("E","1.5E2",0.20,126.45136))
    push!(psvTable, ("F","1.5F3",0.31,198.06412))
    push!(psvTable, ("G","1.5G3",0.50,324.51548))
    push!(psvTable, ("H","2H3",0.79,506.4506))
    push!(psvTable, ("J","3J4",1.29,830.32092))
    push!(psvTable, ("K","3K4",1.84,1185.80408))
    push!(psvTable, ("L","4L6",2.85,1840.64148))
    push!(psvTable, ("M","4M6",3.60,2322.576))
    push!(psvTable, ("N","4N6",4.34,2799.9944))
    push!(psvTable, ("P","4P6",6.38,4116.1208))
    push!(psvTable, ("Q","6Q8",11.05,7129.018))
    push!(psvTable, ("R","6R10",16.00,10322.56))
    push!(psvTable, ("T","8T10",26.00,16774.16))

    maxIndex = size(psvTable.Designation)[1]
    i = 1
    while ((psvTable.areaMM2[i] < area) & ( i < maxIndex))
        i += 1
    end
    
    return (psvTable[i,:Designation])
end

function thermExpansionRate(heat, alpha, heatCap)
    # liquid thermal expansion sizing
    # heat kJ/s,
    # alpha thermal expansion coefficient, 1/K
    # heatCap kJ/kg.K
    # return mass flow(kg/hr) = 3600 * alpha (1/K) * heat (kJ/s) / (heatCap (kJ/kg.K))
    return (3600.0 * alpha * heat / heatCap)
end

function liquidVaporizeReliefRate(heat,latent)
    # heat kJ/s
    # latent is in kJ/kg
    # flow rate kg/h = (heat/1000) * (1/latent) * (1/3600)

    flowRate = heat * (1/latent)* 3600
    return (flowRate)
end

function poolFireReliefRate(wettedAreaM2,latent,prompt)
    # wettedAreaM2 is in m2
    # latent is in kJ/kg
    # prompt fire response is either "prompt" or something else
    # heat Watts = C1 F A^0.82
    # flow rate kg/h = (heat/1000) * (1/latent) * (1/3600)

    C1 = 70900.0   # prompt fire fighting DOES NOT exist

    if (prompt == "prompt")
        C1 = 43200.0   # prompt fire fighting exists
    end
    
    F = 1.0 # no credits
    Q = C1*F*wettedAreaM2^0.82
    flowRate = (Q/1000) * (1/latent)* 3600
    return (flowRate)
end

function PSVsteamRate(areaMM2, Pkpa, State)
    # Napier equation for steam flow through PSV
    # Good for choked flow of steam through orifice
    # Pkpa is in kPaa
    # our function will return the flow rate in kg/hr
    # flowRate = Area * Flux
    psvFlux = PSVsteamFlux(Pkpa, State) # get the flux through the PSV, kg/hr.mm2
    Wkg = areaMM2*psvFlux
    return (Wkg)
end

function PSVsteamSize(Wkg, Pkpa, State)
    # Napier equation for steam flow through PSV
    # Good for choked flow of steam through orifice
    # our function will return the area in MM2
    # flowRate = Area * Flux
    psvFlux = PSVsteamFlux(Pkpa, State) # get the flux through the PSV, kg/hr.mm2
    Amm2 = Wkg/psvFlux
    return (Amm2)
end


# steam flux needs a field for State: Saturated or a temperature in C
function PSVsteamFlux(Pkpa, State)
    # Napier equation for steam flow through PSV
    # mass flow = area*flux
    # flux will be in kg/(mm2.hr)
    # Good for choked flow of steam through orifice
    # Wlbhr = Ain2 Cnapier Kd Ppsi Ksh Kb Kn
    # Wkghr = (1/2.205)*(1/25.4^2)*Amm2 Cnapier Kd Ppsi Ksh Kb Kn
    # flux = (1/2.205)*(1/25.4^2) Cnapier Kd Ppsi Ksh Kb Kn
    # our function will return the flow in kg/h
    # is the napier constant 51.45 or something else?
    # must agree on the value for Kd, 0.975 or 0.938
    # discharge coefficient is usually 0.975 for steam and vapour but
    # could be different depending on manufacturer
    # flowRate = Area * Flux
    
    Cnapier = 51.45
    Kd = 0.975
    Ksh = getKsh(Pkpa, State)
    Kb = 1.0
    Kn = 1.0     # low pressure and high pressure correction
    if (Pkpa > 10300.0)
        Kn = (2.7644*Pkpa/100.0 - 1000.0)/(3.3242*Pkpa/100.0 - 1061.0)
    end
    
    Ppsi = Pkpa * (14.503773800721813/100)
    return (Cnapier*Kd*Ppsi*Ksh*Kb*Kn/(2.205*25.4^2))

end

function PSVvaporFlux(P, Tcelcius, MW, k, Z)
    # vapour equation for PSV rating
    # P in kPaa
    # Tcelcius in celcius,
    # k is ratio of specific heats
    # Z compressibility factor
    # MW mole weight
    # Patm is atm pressure in kPaa
    # return value is flux in kg/mm2.hr
#    coeffSI = 13160.0; choop this out, not in the new version
    T = Tcelcius + 273.15;
    Kd = 0.975 # this could be reduced
    Kb = 1.0 # do not consider backpressure derating
    Kc = 1.0 # no derating for rupture disc
    C = apiC(k)
    rootTerm = sqrt(T*Z/MW);
    coeffTerm = 1.0 / (C * Kd * P * Kb * Kc)
    #    areaMM2 = W*coeffTerm*rootTerm;
    # W = areaMM2 * flux
    # flux = 1.0/(CoeffTerm * rootTerm)

    return (1.0/(coeffTerm*rootTerm))
end

function apiC(k)
    # SI form
    C = 0.03948*sqrt(k*(2.0/(k+1.0))^((k+1.0)/(k-1.0)))   # C coefficient API 520A fig 32
    return (C)
end


function PSVvaporRate(areaMM2, P, Tcelcius, MW, k, Z)
    # vapour equation for PSV rating
    # areaMM2 orifice area mm2
    # W flow rate in kg/h
    # P in kPaa
    # Tcelcius in celcius,
    # k is ratio of specific heats
    # Z compressibility factor
    # MW mole weight
    # Patm is atm pressure in kPaa
    # return value is rated flow rate kg/h
    # W = areaMM2 * flux
    # flux = 1.0/(CoeffTerm * rootTerm)
    psvFlux = PSVvaporFlux(P, Tcelcius, MW, k, Z) # get the flux through the PSV, kg/hr.mm2

    return (areaMM2 * psvFlux)
end

function PSVvaporSize(W, P, Tcelcius, MW, k, Z)
    # vapour equation for PSV sizing
    # W flow rate in kg/h
    # P in kPaa
    # Tcelcius in celcius,
    # k is ratio of specific heats
    # Z compressibility factor
    # MW mole weight
    # Patm is atm pressure in kPaa
    # return value is required area mm2
    # W = areaMM2 * flux
    # flux = 1.0/(CoeffTerm * rootTerm)
    psvFlux = PSVvaporFlux(P, Tcelcius, MW, k, Z) # get the flux through the PSV, kg/hr.mm2

    return (W / psvFlux)
end

function PSVliquidSize(W, P, Pback, d, mu)
    # liquid equation for PSV sizing, PSV requiring capacity certification
    # W flow rate in kg/h
    # P in kPag
    # Pback in kPag
    # d density kg/m3
    # mu viscosity cP or mPa.s
    # Patm is atm pressure in kPaa
    # return value is required area mm2
    coeffSI = 11.78;
    Q = 1000.0*W/(d*60.0);  # litres/minute
    G = d/1000.0; # specific gravity
    P1 = P;
    P2 = Pback;
    
    Kd = 0.65 # This is for a PSV. If we look for a rupture disk, Kd = 0.62
    Kw = 1.0 # do not consider backpressure derating
    Kc = 1.0 # no derating for rupture disc
    
    Kv = 1.0 # first pass at viscous correction, this will be iterated
    i = 1
    areaMM2 = 0.0;
    maxIter = 10;
    while ( i < maxIter)      # iterate 10 times, we are not going to check for convergence
        rootTerm = sqrt(G/(P1-P2));
        coeffTerm = coeffSI * Q / (Kd * Kw * Kc * Kv)
        areaMM2 = coeffTerm*rootTerm;
#        diam = sqrt(4.0*areaMM2/pi)/1000;
        R = Q * 18800 * G / (mu*sqrt(areaMM2)); # Reynolds number
        Kv = 1.0 / (0.9935 + 2.878/sqrt(R) + 342.75/(R^1.5))

        i += 1
    end
    return (areaMM2)
end



function PSVliquidRate(areaMM2, P, Pback, d, mu)
    # liquid equation for PSV sizing
    # W flow rate in kg/h
    # areaMM2 PSV area
    # P in kPag
    # Pback in kPag
    # d density kg/m3
    # mu viscosity cP or mPa.s
    # Patm is atm pressure in kPaa
    # return value is flow rate kg/h
    # I have not edited this function yet
    # from API 520
    # A = (11.78 Q) / (Kd Kw Kc Kv) * sqrt(G/DP)
    # Q = (A Kd Kw Kc Kv/11.78) sqqrt(DP/G)
    coeffSI = 11.78;
#    Q = 1000*W/(d*60.0);  # litres/minute
#    A = areaMM2 / (25.4^2)
    G = d/1000; # specific gravity
    P1 = P;
    P2 = Pback;

    Kd = 0.65 # typical value
    Kw = 1.0 # do not consider backpressure derating
    Kc = 1.0 # no derating for rupture disc

    Kv = 1.0 # first pass at viscous correction, this will be iterated
    i = 1
    Q = 0.0 # initialize the variable
    W = 0.0;  # initialize the variable
    maxIter = 10;
    while ( i < maxIter)                  # iterate 10 times, we are not going to check for convergence
        # areamm2 = coeff * root
        rootTerm = sqrt((P1-P2)/G);
        coeffTerm = (Kd * Kw * Kc * Kv) / coeffSI
        Q = areaMM2*coeffTerm*rootTerm;   # litres per minute
#        diam = sqrt(4.0*areaMM2/pi)/1000.0;
        R = Q * 18800.0 * G / (mu*sqrt(areaMM2));  # Reynolds number
        Kv = 1.0 / (0.9935 + 2.878/sqrt(R) + 342.75/(R^1.5))
        
        i += 1
    end
    return (Q*60*d/1000.0)  # mass flow rate kg/h
end

end # module
