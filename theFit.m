classdef theFit < audioPlugin
    %-----------------------------------------------------------------------
    % Private Properties - Used for internal storage
    %-----------------------------------------------------------------------
    properties (Access = private)
        %Modeling characteristics
        asym_n = 1;
        asym_p = 1;
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
        Ts = 1/Fs;
        localMax = zeros(10, 2);
        Is_neg = Is / (asym_n * eta * Vt);
        Is_pos = Is / (asym_p * eta * Vt);
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
                'Mapping',{'lin',0,4})...
                );
        end
    
    methods
        %-------------------------------------------------------------------
        % Main processing function
        %-------------------------------------------------------------------
        function y = process(self, x)
            [N, channels] = size(x);
            %Retrieve private variable to use locally in loop
            Is_neg = self.Is_neg;
            Is_pos = self.Is_pos;
            
        end
        
        %-------------------------------------------------------------------
        % Set Method
        %-------------------------------------------------------------------
        function set.Cutoff(self, val)
            
        end
        
        %-------------------------------------------------------------------
        % Reset Method
        %-------------------------------------------------------------------
        function reset(self)
            
            
        end
        %-------------------------------------------------------------------
        % Set Methods for Harmonics
        %-------------------------------------------------------------------
        function set.second(self, in)
            setDiodes('2');
        end
        function set.third(self, in)
            setDiodes('3');
        end
        function set.fourth(self, in)
            setDiodes('4');
        end
        function set.fifth(self, in)
            setDiodes('5');
        end
        function set.sixth(self, in)
            setDiodes('6');
        end
        function set.seventh(self, in)
            setDiodes('7');
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
        % Interpret user settings
        %-------------------------------------------------------------------
        function [] = setDiodes(self, numDiode)
            if(in ~= self.NDiodes)
                return
            else
                %If there are more neg diodes then pos diodes
                %
                if(in > self.PDiodes)
                    
                %If there are more pos diodes then neg diodes   
                elseif(in < self.PDiodes)
                        
                end
            end
            
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
                    self.harms(j) = self.localMax(i, 1);
                    j = j + 1;
                end
                if currMaxFreq == self.localMax(i, 1)
                    self.harms(1) = currMaxFreq;
                    found = true;
                    j = 2;
                end
            end
        end
    end
end