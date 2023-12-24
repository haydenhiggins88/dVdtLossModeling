
%does positive / negative Load Current matter? Do I need to account for
%both Th and Tl in my calculations??

%Why are Rint and Cgd different for turn-on and turn-off (see 3)experimental validation)? which values do I
%use?

Vgon = 15;
Vgoff = -6;
Vth = 6;
global Rgon;                %30.1 16.2 16.2      15 8.6 5.2   %VARIABLES UNDER TEST:
global Rgoff;                 %24.3 16.2 11.0     15 10 6.8
global Cm;              %20 50    pF
Cgdon = 24.5*10^-12;          %for turn-on, 45.7 for turn off??   pF
Cgdoff = 45.7*10^-12;
Cds = 230*10^-12;           %Coss for C3M0016120K       pF
Rinton = 2.2;                 %for turn on, 2.5 for turn off  ??
Rintoff = 2.5;
TAUon = 143*10^-12;         %ps
TAUoff = 151*10^-12;        %ps
Vm = 8.7;                   %Miller voltage ?
Ceff = 666*10^-12;          %pF, = 2CdQ,oss
Vdc = 800;
Qtot = 384*10^-9;           %nC

%dVdtZCS = 7.5*10^9;        %LUT w/ respect to miller capacitance? see figure 5

Fsw = 100*10^3 ;            %Switching Frequency
Fsin = 50 ;                 %Sinusoid frequency ; 50 Hz = 3000 RPM
Nsw = (Fsw / Fsin)/4;           %Number of switching events per quarter cycle, 2000 in this case

Isin = zeros(Nsw,1);        %blank LUT array
for i = 1 : Nsw
    i;
    Isin(i) = round (45*sin(2*pi*i / (4*Nsw) ) );
end
plot(Isin)                  %LUT for A=45 sinusoid, resolution = 2000
Esw = zeros(46,1);          %LUT for loss of a switching event at particular Isw
Elosses1 = zeros(Nsw,1);    %table to store losses from this particular simulation
Elosses2 = zeros(Nsw,1);
Elosses3 = zeros(Nsw,1);
Power = zeros(3,1)          %


% -----     SIM WITH ONE SETTING     ------ %

for n = 1:Nsw
    if n == 1
        setParam(16.2, 11, 0)
        Rgon 
        Rgoff 
        Cm 
    end
    if n == 1 %recalculates constants when Rg and Cm get updated (different "windows" throughout the Isw sinusoid)
        dVdtoff = (Vm - Vgoff) / (Rgoff*Cm + Rgoff*Cgdoff + Rintoff*Cgdoff + TAUoff);
        dVdton = (Vgon - Vth) / (Rgon*Cm + Rgon*Cgdon + Rinton*Cgdon + TAUon);
        Ik = (1 + Cds / Cgdoff)*(Vth + Vgoff) / Rgoff;      %Ik = Ceff*dVdtZCS
        koff = ( (1/2) * Vdc^2) / dVdtoff;
        kon = ( (1.35/2) * Vdc^2) / dVdton;

        for i = 1 : 46 %matlab is stupid and doesn't let you index from 0, so that's why there's a lot of off-by-one corrections
            Isw = i - 1;
            if Isw < Ik
              Isw;
              Esw(n) = Vdc*Qtot + kon*Iswitch;
            end
            if Isw > Ik
                Isw;
                Esw(i) = Vdc*Qtot + kon*Isw + koff*(Isw - Ik);
            end
        end
    end

    Elosses1(n) = Esw(abs(Isin(n)) + 1);   %the indices of Esw are +1, e.g. to find Esw(0) you must go to index 1
end

% -----      SIM WITH TWO SETTINGS ------ %

for n = 1:Nsw
    if n == 1
        setParam(16.2, 11, 0)
        Rgon
        Rgoff
        Cm
    end
    if n == Nsw/2
        setParam(10, 5, 0)
        Rgon
        Rgoff
        Cm
    end
    if n == 1 || n == Nsw/2   %recalculates constants when Rg and Cm get updated (different "windows" throughout the Isw sinusoid)
        dVdtoff = (Vm - Vgoff) / (Rgoff*Cm + Rgoff*Cgdoff + Rintoff*Cgdoff + TAUoff);
        dVdton = (Vgon - Vth) / (Rgon*Cm + Rgon*Cgdon + Rinton*Cgdon + TAUon);
        Ik = (1 + Cds / Cgdoff)*(Vth + Vgoff) / Rgoff;      %Ik = Ceff*dVdtZCS
        koff = ( (1/2) * Vdc^2) / dVdtoff;
        kon = ( (1.35/2) * Vdc^2) / dVdton;

        for i = 1 : 46 %matlab is stupid and doesn't let you index from 0, so that's why there's a lot of off-by-one corrections
            Isw = i - 1;
            if Isw < Ik
              Isw;
              Esw(n) = Vdc*Qtot + kon*Iswitch;
            end
            if Isw > Ik
                Isw;
                Esw(i) = Vdc*Qtot + kon*Isw + koff*(Isw - Ik);
            end
        end
    end

    Elosses2(n) = Esw(abs(Isin(n)) + 1);   %the indices of Esw are +1, e.g. to find Esw(0) you must go to index 1
end

% -----      SIM WITH THREE SETTINGS ------ %

for n = 1:Nsw
    if n == 1
        setParam(16.2, 11, 0)
        Rgon
        Rgoff
        Cm
    end
    if n == round(Nsw/3)
        setParam(10, 5, 0)
    end
    if n == round((2*Nsw/3))
        setParam(5, 3, 0)
        Rgon
        Rgoff
        Cm
    end
    if n == 1 || n == round(Nsw/3)  || n == round((2*Nsw/3));   %recalculates constants when Rg and Cm get updated (different "windows" throughout the Isw sinusoid)
        dVdtoff = (Vm - Vgoff) / (Rgoff*Cm + Rgoff*Cgdoff + Rintoff*Cgdoff + TAUoff);
        dVdton = (Vgon - Vth) / (Rgon*Cm + Rgon*Cgdon + Rinton*Cgdon + TAUon);
        Ik = (1 + Cds / Cgdoff)*(Vth + Vgoff) / Rgoff;      %Ik = Ceff*dVdtZCS
        koff = ( (1/2) * Vdc^2) / dVdtoff;
        kon = ( (1.35/2) * Vdc^2) / dVdton;

        for i = 1 : 46 %matlab is stupid and doesn't let you index from 0, so that's why there's a lot of off-by-one corrections
            Isw = i - 1;
            if Isw < Ik
              Isw;
              Esw(n) = Vdc*Qtot + kon*Iswitch;
            end
            if Isw > Ik
                Isw;
                Esw(i) = Vdc*Qtot + kon*Isw + koff*(Isw - Ik);
            end
        end
    end

    Elosses3(n) = Esw(abs(Isin(n)) + 1);   %the indices of Esw are +1, e.g. to find Esw(0) you must go to index 1
end

%Elosses1
%Elosses3
Ecycle(1) = 4*sum(Elosses1)                   %simulation solves for the losses from a quarter of a sine period, so *4 --> full period
Ecycle(2) = 4*sum(Elosses2)
Ecycle(3) = 4*sum(Elosses3)

Power = Ecycle*Fsin


% -----     FUNCTIONS     ------ %

function [] = setParam(newRgon, newRgoff, newCm)
    global Rgon;
    global Rgon;
    global Rgoff

    Rgon = newRgon;
    Rgoff = newRgoff;
    Cm = newCm;
end