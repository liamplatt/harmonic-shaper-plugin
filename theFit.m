classdef theFit < audioPlugin
    %-----------------------------------------------------------------------
    % Private Properties - Used for internal storage
    %-----------------------------------------------------------------------
    properties (Access = private)
        %Modeling characteristics
        Ra = 1000;
        R6 = 51e3;
        R4 = 4.7e3;
        R5 = 10e3;
        C4 = 90e-12;
        C3 = .047e-6;
        Voff = 4.5;
        prec = 1e-6;
        maxIter = 50;
        Vt = 0.026;
        eta = 1;
        Is = 10^-12;
        localMax = zeros(10, 2);
        Is_neg; 
        Is_pos;
        %Set initial state variables
        b0 = 0.5;
        b1 = 0.5;
        a1 = 0.1;
        x = 0;
        xhold = zeros(1, 1);
        yhold = zeros(1, 1);
        g = 0.5;
        pBins = [20,39, 40,79, 80,159, 160,319, 320,639, 640 1279 1280 2559 2560 5119  5120 10239 10240 20000];
        plocalMax = zeros(10, 2);
        pHarms = zeros(7,1);
        pPDiodes = 1;
        pNDiodes = 1;
    end
    %-----------------------------------------------------------------------
    % Public Properties - End user interacts with these
    %-----------------------------------------------------------------------
    properties
        second = 1;
        third = 1;
        fourth = 1;
        fifth = 1;
        sixth = 1;
        seventh = 1;
        drive = 1;
    end
    %-----------------------------------------------------------------------
    % Constant Properties - Used to define plugin interface
    %-----------------------------------------------------------------------
    properties (Constant)
        PluginInterface = audioPluginInterface( ...
            audioPluginParameter('second', ...
                'Label','2nd Harm', ...
                'Mapping',{'lin',0,4}),...
            audioPluginParameter('third', ...
                'Label','3rd Harm', ...
                'Mapping',{'lin',0,4}),...
            audioPluginParameter('fourth', ...
                'Label','4th Harm', ...
                'Mapping',{'lin',0,4}),...
            audioPluginParameter('fifth', ...
                'Label','5th Harm', ...
                'Mapping',{'lin',0,4}),...
            audioPluginParameter('sixth', ...
                'Label','6th Harm', ...
                'Mapping',{'lin',0,4}),...    
            audioPluginParameter('seventh', ...
                'Label','7th Harm', ...
                'Mapping',{'lin',0,4}),...
            audioPluginParameter('drive', ...
                'Label','gain', ...
                'Mapping',{'log',0.000001,6})...
                );
        end
    
    methods
        %-------------------------------------------------------------------
        % Main processing function
        %-------------------------------------------------------------------
        function Vout = process(self, Vin)
            [N, channels] = size(x);
            Vout = zeros(N, channels);
            %Retrieve private variable to use locally in loop
            Is_neg = self.Is_neg;
            Is_pos = self.Is_pos;
            eta = self.eta;
            Vt = self.Vt;
            x = self.x;
            R = self.R;
            Is = self.Is;
            b0 = self.b0;
            b1 = self.b1;
            a1 = self.a1; 
            maxIter = self.maxIter;
            prec = self.prec;
            for c = 1:channels
                for n = 1:N
                    i = 1;          %positive diode                                                             %negative diode
                    fx = x/R + Is * (exp(x/(asym_p * eta*Vt))-1) - Vin(n,channel)/R + -Vout(n,channel)/R - Is * (exp(-x/(asym_n*eta*Vt))-1);
                    %Solve for voltage at each sample
                    while(i < maxIter && (abs(fx)> prec))
                        %Take dirivitive to and solve using newton raphson method
                        den = 1/R + Is_pos * (exp(x/(asym_p * eta*Vt)))+ Is_neg * (exp(-x/(asym_n*eta*Vt)));
                        x = x - (fx/den);
                        i = i + 1;
                        fx = x/R + Is * (exp(x/(asym_p * eta*Vt))-1) - Vin(n,channel)/R + -Vout(n,channel)/R - Is * (exp(-x/(asym_n*eta*Vt))-1);
                    end
                    Vout(n,channel) =  b0 * x + b1 * xhold - a1 * yhold;
                    %update state variables
                    yhold = Vout(n, channel);
                    xhold = x;
                end
            end
        end
        %-------------------------------------------------------------------
        % Reset Method
        %-------------------------------------------------------------------
        function reset(self)
            self.plocalMax = zeros(10, 2);
            self.pHarms = zeros(7,1);
            self.pPDiodes = 1;
            self.pNDiodes = 1;
        end
        %-------------------------------------------------------------------
        % Set Methods for Harmonics
        %-------------------------------------------------------------------
        function set.second(self, in)
            self.setDiodes(2);
        end
        function set.third(self, in)
            self.setDiodes(3);
        end
        function set.fourth(self, in)
            self.setDiodes(4);
        end
        function set.fifth(self, in)
            self.setDiodes(5);
        end
        function set.sixth(self, in)
            self.setDiodes(6);
        end
        function set.seventh(self, in)
            self.setDiodes(7);
        end
    end
    methods (Access = private)
        %-------------------------------------------------------------------
        % Calculate Filter Coefficients
        %-------------------------------------------------------------------
        function [B,A] = notchCoeffs(~,fc,fs)
            [B, A] = iirnotch(fc/fs, .01, 3); 
        end
        %-------------------------------------------------------------------
        % Set function for private props pPDiodes pNDiodes
        %-------------------------------------------------------------------
        function [] = setAssym(self, neg, pos)
           self.pPDiodes = pos;
           self.pNDiodes = neg;
           reset(self);
        end
        %-------------------------------------------------------------------
        % Helper function for swapping properties
        %-------------------------------------------------------------------
        function [] = swapNumDiodes(self)
            temp = self.pNDiodes;
            self.pNDiodes = self.pPDiodes;
            self.pPDiodes = temp;
        end
        %-------------------------------------------------------------------
        % Pos=Even Neg=Odd
        %-------------------------------------------------------------------
        function [] = setTwo(self)
            %Take into count current assymp and assymneg
            if self.pPDiodes < self.pNDiodes
                swapNumDiodes(self);
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setThree(self)
            if self.pPDiodes < self.pNDiodes
                swapNumDiodes(self);
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setFour(self)
            if self.pPDiodes < self.pNDiodes
                swapNumDiodes(self);
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setFive(self)
            if self.pPDiodes < self.pNDiodes
                swapNumDiodes(self);
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setSix(self)
            if self.pPDiodes < self.pNDiodes
                swapNumDiodes(self);
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setSeven(self)
            if self.pPDiodes < self.pNDiodes
                swapNumDiodes(self);
            end
        end
        %-------------------------------------------------------------------
        % Set Diodes, helper function, calls each of the harms function
        %-------------------------------------------------------------------
        function [] = setDiodes(self, numHarm)
            switch numHarm
                case 2
                    setTwo(self);
                case 3
                    setThree(self);
                case 4
                    setFour(self);
                case 5
                    setFive(self);
                case 6
                    setSix(self);
                case 7
                    setSeven(self);
                otherwise
                    return;
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setHold(self, xhold, yhold)
            self.xhold = xhold;
            self.yhold = yhold;
        end
        %-------------------------------------------------------------------
        % Calculate local fft peaks
        %-------------------------------------------------------------------
        function [] = calcHarms(self, in)
            xdft = fft(Vout);
            xdft = xdft(1:floor(N/2+1));
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            psdx = 10 * log(psdx);
            
            for i = 1:2:20
                currentBin = self.bins(i:i+1);
                currMaxAmp = -inf;
                currMaxFreq = 1;
                for j = currentBin(1):currentBin(2)
                    psdx(j);
                    if psdx(j) > currMaxAmp
                        currMaxAmp = psdx(j);
                        currMaxFreq = j;
                    end
                end
                self.localMax(k, 1) = currMaxFreq;
                self.localMax(k, 2) = currMaxAmp;
                k = k + 1;
            end
            %once we have local maximums, then we go through and find maximum there
            currMaxAmp = -inf;
            found = false;
            for i = 1:10
                if self.localMax(i, 2) > currMaxAmp
                    currMaxAmp = self.localMax(i, 2);
                    currMaxFreq = self.localMax(i, 1);
                end
            end
            for i = 1:10
                if found
                    self.pHarms(j) = self.localMax(i, 1);
                    j = j + 1;
                end
                if currMaxFreq == self.localMax(i, 1)
                    self.pHarms(1) = currMaxFreq;
                    found = true;
                    j = 2;
                end
            end
        end
    end
end