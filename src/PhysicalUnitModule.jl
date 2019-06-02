"""
    PhysicalUnitModule

Module  to handle use of physical units in definitions of input data.
"""
module PhysicalUnitModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
# using Unicode


"""
    phun(str::String; system_of_units = :SI, base_time_units = :SEC)

Evaluate an expression in physical units.

# Example
```
pu = ustring -> phun(ustring; system_of_units = :SIMM)
E1s = 130.0*pu("GPa")
```
yields
1.3e+5 (in mega Pascal)
whereas
```
130.0*phun("GPa"; system_of_units = :SI)
```
yields
1.3e+11 (in Pascal)
"""
function phun(str::String; system_of_units = :SI, base_time_units = :SEC)::FFlt

    @assert system_of_units in [:SI :SIM :US :IMPERIAL :CGS :SIMM]
    @assert base_time_units in [:SEC :MIN :HR :DY :YR :WK]

    """
        _physunitdict(system_of_units, base_time_units)
    Create a dictionary with physical unit conversion factors.
    Inputs:
    --`system_of_units`
    if system_of_units  ==  :US
       basic assumed units are American Engineering:
       LENGTH = FT, TIME = SEC, MASS = SLUG TEMPERATURE = RAN FORCE = LB
    elseif system_of_units  ==  :CGS
       basic assumed units are Centimeter,Gram,Second:
       LENGTH = CM, TIME = SEC, MASS = GM TEMPERATURE = K FORCE = DYNE
    elseif system_of_units  ==  :IMPERIAL
       basic assumed units are Imperial:
       LENGTH = FT, TIME = SEC, MASS = SLUG TEMPERATURE = RAN FORCE = LB
    otherwise,
       basic assumed units are :SIM (equivalent to :SI, default):
       LENGTH = M , TIME = SEC, MASS = KG   TEMPERATURE = K   FORCE = N
    --`base_time_units` defaults to :SEC
    """
    function _physunitdict(system_of_units, base_time_units)

        d =  Dict{String, FFlt}();

        uSEC  =  1.0;
        if base_time_units  ==  :MIN
            uSEC  =  1.0/60;
        elseif base_time_units  ==  :HR
            uSEC  =  1.0/(60*60);
        elseif base_time_units  ==  :DY
            uSEC  =  1.0/(60*60*24);
        elseif base_time_units  ==  :YR
            uSEC  =  1.0/(60*60*24*365);
        elseif base_time_units  ==  :WK
            uSEC  =  1.0/(60*60*24*7);
        end
        d["SEC"] = uSEC;

        # time conversions
        uMIN  =  uSEC*60; d["MIN"] = uMIN;
        uHR  =  uSEC*60*60; d["HR"] = uHR;
        uDAY  =  uSEC*60*60*24; d["DAY"] = uDAY;
        uWK  =  uSEC*60*60*24*7; d["WK"] = uWK;
        uMONTH  =  uSEC*60*60*24*365/12; d["MONTH"] = uMONTH;
        uYR  =  uSEC*60*60*24*365; d["YR"] = uYR;

        # base units for charge, time, angle (all systems)
        uCOUL = 1.0;uRAD = 1.0; d["COUL"] = uCOUL; d["RAD"] = uRAD;

        # Base units for length, mass, temperature (different for each system)
        # :SI or :SIM (SI units using meters)
        uM = 1.0; uKG = 1.0; uK = 1.0;						  # SI base units
        uCM = uM/100.0; uGM = uKG/1000.0;               # CGS base units
        uFT = uM/(3*(1/0.9144)); uSLUG = uKG*14.593903591998501; uRAN = uK/1.8; # US base units
        if system_of_units  ==  :US
            uFT = 1.0; uSLUG = 1.0; uRAN = 1.0;                     # US base units
            uM = uFT*3*(1/0.9144); uKG =  uSLUG/14.593903591998501; uK = 1.8*uRAN; # SI base units
            uCM = uM/100.0;uGM =  uKG/1000.0;                     # CGS base units
        elseif system_of_units  ==  :IMPERIAL
            uFT = 1.0; uSLUG = 1.0; uRAN = 1.0;                     # US base units
            uM = uFT*3*(1/0.9144); uKG =  uSLUG/14.593903591998501; uK = 1.8*uRAN; # SI base units
            uCM = uM/100.0;uGM =  uKG/1000.0;                     # CGS base units
        elseif system_of_units  ==  :CGS
            uCM = 1.0; uGM = 1.0; uK = 1.0;						  # CGS base units
            uM = 100.0*uCM;uKG = 1000.0*uGM; 					  # SI base units
            uFT = uM/(3*(1/0.9144)); uSLUG = uKG*14.593903591998501; uRAN = uK/1.8; # US base units
        elseif system_of_units  ==  :SIMM #SI(mm)
            uM = 1000.0; uKG = 1.0e-3; uK = 1.0;						  # SI base units
            uCM = uM/100.0; uGM = uKG/1000.0;             # CGS base units
            uFT = uM/(3*(1/0.9144)); uSLUG = uKG*14.593903591998501; uRAN = uK/1.8; # US base units
        end
        d["M"]   = uM;
        d["KG"] = uKG;
        d["K"] = uK;
        d["CM"] = uCM;
        d["GM"] = uGM;
        d["FT"] = uFT;
        d["SLUG"] = uSLUG;
        d["RAN"] = uRAN;


        # prefixes
        uNANO = 10.0^(-9); d["NANO"] = uNANO;
        uMICRO = 10.0^(-6); d["MICRO"] = uMICRO;
        uMILLI = 10.0^(-3); d["MILLI"] = uMILLI;
        uKILO = 10.0^3; d["KILO"] = uKILO;
        uMEGA = 10.0^6; d["MEGA"] = uMEGA;
        uGIGA = 10.0^9; d["GIGA"] = uGIGA;
        uTERA = 10.0^12; d["TERA"] = uTERA;

        # angle conversions (degrees of angle)
        uDEG = pi*uRAD/180.0; d["DEG"] = uDEG;
        uREV = 2.0*pi;# Measure of angles in terms of revolutions
        d["REV"] = uREV;

        # length conversions
        uIN = uFT/12; d["IN"] = uIN;
        uNMI = 1852*uM; d["NMI"] = uNMI;
        uCM = uM/100.0; d["CM"] = uCM;
        uYD = 3.0*uFT; d["YD"] = uYD;
        uMM = uM/1000.0; d["MM"] = uMM;
        uMILE = 5280*uFT; d["MILE"] = uMILE;

        # acceleration conversions
        uG = 32.17405*uFT/uSEC^2; d["G"] = uG;

        # force conversions
        uLBF = uSLUG*uFT/uSEC^2; d["LBF"] = uLBF;
        uNT = uKG*uM/(uSEC^2); d["NT"] = uNT;  #newton
        uOZ = uLBF/16.0; d["OZ"] = uOZ;

        # mass conversions
        uLBM = uLBF/(uG); d["LBM"] = uLBM;

        # speed conversions
        uKT = uNMI/uHR; d["KT"] = uKT;
        uMPH = uMILE/uHR; d["MPH"] = uMPH;

        # pressures
        uPSI = uLBF/(uIN^2); d["PSI"] = uPSI;
        uPA =  uNT/(uM^2); d["PA"] = uPA;
        uBAR = 10^5*uPA; d["BAR"] = uBAR;
        uATM = 14.6959*uPSI; d["ATM"] = uATM;
        uTORR = 133.322*uPA; d["TORR"] = uTORR;
        ummHG = uTORR; d["MMHG"] = ummHG;
        uBA = 0.1*uPA;  d["BA"] = uBA;#cgs unit of pressure, Barye

        # work
        uJ =  uNT*uM; d["J"] = uJ;
        uCAL  =  uJ * 4.1868; d["CAL"] = uCAL;
        uMEV  =  uJ / 6.242e12; d["MEV"] = uMEV;
        uERG  =  uJ * 1E-7; d["ERG"] = uERG;
        uBTU  =  1054*uJ; d["BTU"] = uBTU;

        # power
        uW = uJ/uSEC; d["W"] = uW;
        uMW  =  1e6*uW; d["MW"] = uMW;
        uHP = 745.7*uW; d["HP"] = uHP;

        # electricity
        uA  =  uCOUL/uSEC; d["A"] = uA;
        uV  =  uJ/uCOUL; d["V"] = uV;
        uOHM  =  uA/uV; d["OHM"] = uOHM;

        # frequency
        uHZ  =  1.0/uSEC; d["HZ"] = uHZ;
        uRPS  =  uHZ; d["RPS"] = uRPS;# Revolutions per second
        uRPM  =  1.0/uMIN; d["RPM"] = uRPM;# Revolutions per minute

        # volume
        uL  =  1000.0*uCM^3; d["L"] = uL;;

        if system_of_units  ==  :IMPERIAL
            uGAL  =  4.54609*uL;# Imperial gallon
        else
            uGAL  =  3.785411784*uL;# US gallon
        end
        d["GAL"] = uGAL;
        uBBL  =  159*uL;# standard barrel for measuring volume of oil http://en.wikipedia.org/wiki/Barrel_(unit)
        d["BBL"] = uBBL;

        # flow rate
        uGPM  =  uGAL/uMIN;# gallons per minute
        d["GPM"] = uGPM;

        # Temperature, degrees Fahrenheit
        uF = uRAN;
        d["F"] = uF;

        # Aliases
        # Alias for a second
        d["S"] = uSEC;
        # The symbol for the Newton
        d["N"] = uNT;
        # The symbol for the mega Pascal
        d["MPA"] = uMEGA*uPA;
        # The symbol for the giga Pascal
        d["GPA"] = uGIGA*uPA;
        # The symbol for the thousands of pounds force
        d["KIPS"] = uKILO*uLBF;
        # The symbol for the thousands of pounds per square inch
        d["KSI"] = uKILO*uPSI;
        # The symbol for the centiPoise
        d["CP"] = 1.0/1000.0*uPA*uSEC;

        return d
    end

    d  =  _physunitdict(system_of_units, base_time_units);

    function replacesymbols(d::Dict{String, FFlt}, str::String)
        str = uppercase(str);
        outstr = "";
        i = 1;
        while i <= length(str)
            if isuppercase(str[i])
                k = i+1;
                while (k <= length(str)) && (isuppercase(str[k]))
                    k = k+1;
                end
                outstr = outstr * string(d[str[i:k-1]]);
                i = k-1;
            else
                outstr = outstr * string(str[i]);
            end
            i = i+1;
        end
        return outstr
    end

    ostr  =  replacesymbols(d, str);
    val  =  FFlt(Core.eval(Main, Meta.parse(ostr)))
    return val::FFlt
end


end
