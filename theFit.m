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
        Is = 10^-9;
        localMax = zeros(10, 2);
        Is_neg = 0; 
        Is_pos = 0;
        %Set initial state variables
        x = 0;
        xhold = [0, 0; 0, 0];
        yhold = [0, 0; 0, 0];
        g = 0.5;
        eqOn = true;
        pBins = [20,39, 40,79, 80,159, 160,319, 320,639, 640 1279 1280 2559 2560 5119  5120 10239 10240 20000];
        plocalMax = zeros(10, 2);
        pHarms = zeros(7,1);
        pPDiodes = 2;
        pNDiodes = 2;
        B = zeros(1, 3);
        A = zeros(1, 3);
        sampsPerFrame = 1024;
        maxDiodes = 4;
        attenFreq = 3;
        b0 = 0.5;
        b1 = 0.5;
        a1 = 0.1;
        Yhold = zeros(2, 1);
        Xhold = zeros(2, 1);
        State = zeros(2, 1);
        State2 = zeros(2, 2);
        isSetup = false;
        channels = 2;
        builtFilter = true;
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
        PluginInterface = audioPluginInterface(...
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
                'Mapping',{'log',1,10})...
                );
        end
    
    methods
        %-------------------------------------------------------------------
        % Main processing function
        %-------------------------------------------------------------------
        function Vout = process(self, Vin)
            [N, channels] = size(Vin);
            Vout = zeros(N, channels);
            if self.channels ~= channels
               self.setChannels(channels);
            end
            Vin = Vin* self.drive;
            %Retrieve private variable to use locally in loop
            Is_neg = self.Is_neg;
            Is_pos = self.Is_pos;
            eta = self.eta;
            Vt = self.Vt;
            x = self.x;
            R = self.Ra;
            Is = self.Is; 
            maxIter = self.maxIter;
            prec = self.prec;
            asym_p = self.getNumNDiodes()
            asym_n = self.getNumPDiodes()
            isSetup = self.isSetup;
            %If we have set any of the user params, we need to recalc harms
            
            freq = self.pHarms(self.attenFreq)

            xhold = self.Xhold;
            yhold = self.Yhold;
            b0 = self.b0;
            b1 = self.b1;
            a1 = self.a1;
            %apply distortion
            for channel = 1:channels
                for n = 1:N
                    i = 1;
                    fx = x/R + Is * (exp(x/(asym_p * eta*Vt))-1) - Vin(n,channel)/R + -Vout(n,channel)/R - Is * (exp(-x/(asym_n*eta*Vt))-1);
                    %Solve for x
                    while(i < maxIter && (abs(fx)> prec))
                        den = 1/R + Is_pos * (exp(x/(asym_p * eta*Vt)))+ Is_neg * (exp(-x/(asym_n*eta*Vt)));
                        x = x - (fx/den);
                        i = i + 1;
                        fx = x/R + Is * (exp(x/(asym_p * eta*Vt))-1) - Vin(n,channel)/R + -Vout(n,channel)/R - Is * (exp(-x/(asym_n*eta*Vt))-1);
                    end
                    if x >= .5
                       x = .5;
                    elseif x <= -.5
                        x = -.5;
                    end
                    Vout(n,channel) =  b0 * x + b1 * xhold(channel) - a1 * yhold(channel);
                    %update state variables
                    yhold(channel) = Vout(n, channel);
                    xhold(channel) = x;
                end
            end
            if ~isSetup
                self.calcHarms(Vin);
            end
            A = self.A;
            B = self.B;
            State = self.State;
            %Update state vars
            %update yhold
            if self.builtFilter
                [Vout, State] = filter(B, A, Vin, State);
                Size = size(State);
                if Size(2) == 1
                    self.State = State;
                elseif Size(2) == 2
                    self.State2 = State;
                end
            end
        end
        %-------------------------------------------------------------------
        % Reset Method
        %-------------------------------------------------------------------
        function reset(self)
            %self.plocalMax = zeros(10, 2);
            self.pHarms = zeros(7,1);
            self.State = zeros(2, 1);
            self.pPDiodes = 1;
            self.pNDiodes = 1;
            self.isSetup = false;
        end
        %-------------------------------------------------------------------
        % Set Methods for Harmonics
        %-------------------------------------------------------------------
        function set.second(self, in)
            self.second = in;
            self.setDiodes(2, in);
        end
        function set.third(self, in)
            self.third = in;
            self.setDiodes(3, in);
        end
        function set.fourth(self, in)
            self.fourth = in;
            self.setDiodes(4, in);
        end
        function set.fifth(self, in)
            self.fifth = in;
            self.setDiodes(5, in);
        end
        function set.sixth(self, in)
            self.sixth = in;
            self.setDiodes(6, in);
        end
        function set.seventh(self, in)
            self.seventh = in;
            self.setDiodes(7, in);
        end        
    end
    methods (Access = private)
        function setChannels(self, channels)
            self.channels = channels;
            self.Xhold = zeros(channels, 1);
            self.Yhold = zeros(channels, 1);
        end
        %-------------------------------------------------------------------
        % Set Is_neg and Is_pos
        %-------------------------------------------------------------------
        function [Is_neg, Is_pos] = setIs(self)
            self.Is_neg = self.Is / (self.pNDiodes * self.eta * self.Vt);
            self.Is_pos = self.Is / (self.pPDiodes * self.eta * self.Vt);
            Is_pos = self.Is_pos;
            Is_neg = self.Is_neg;
        end
        %-------------------------------------------------------------------
        % Get neg diodes
        %-------------------------------------------------------------------
        function out = getNumNDiodes(self)
            out = self.pNDiodes;
        end
        %-------------------------------------------------------------------
        % Get pos diodes
        %-------------------------------------------------------------------
        function out = getNumPDiodes(self)
            out = self.pPDiodes;
        end
        %-------------------------------------------------------------------
        % Calculate Filter Coefficients
        %-------------------------------------------------------------------
        function [] = buildFiltCoeffs(self, fc)
            Q = 10;
            w0 = (2*pi*fc)/self.getSampleRate();
            alpha = sin(w0)/(2*Q);
            B = [1, -2*cos(w0), 1];
            A = [1+alpha, -2*cos(w0), 1-alpha];
            for i = 1:3
                if isnan(B(i))
                    B(i) = 0;
                elseif isnan(A(i))
                    A(i) = 0;
                end
            end
            self.B = B;
            self.A = A;
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
        function [] = setTwo(self, in)
            %Take into count current assymp and assymneg
            if in > self.pPDiodes
                self.pPDiodes = self.pPDiodes+.1;
                if self.pPDiodes >= 6
                    %Cant increase num diodes anymore
                    self.pPDiodes = 6;
                end
            elseif in == self.pPDiodes
                return;
            else
                self.pPDiodes = self.pPDiodes-.1;
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setThree(self, in)
            if self.pNDiodes >= self.maxDiodes
                self.attenFreq = 3;
            else
                self.pNDiodes = self.pNDiodes+.1;
                if self.pNDiodes >= self.maxDiodes
                    %Cant increase num diodes anymore
                    self.pNDiodes = self.maxDiodes;
                end
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setFour(self, in)
            if in > self.pPDiodes
                self.pPDiodes = self.pPDiodes+.1;
                if self.pPDiodes >= 6
                    %Cant increase num diodes anymore
                    self.pPDiodes = 6;
                end
            elseif in == self.pPDiodes
                return;
            else
                self.pPDiodes = self.pPDiodes-.1;
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setFive(self, in)
            if in > self.pNDiodes
                self.pNDiodes = self.pNDiodes+.1;
                if self.pNDiodes >= 6
                    %Cant increase num diodes anymore
                    self.pNDiodes = 6;
                end
            elseif in == self.pNDiodes
                return;
            else
                self.pNDiodes = self.pNDiodes-.1;
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setSix(self, in)
            if in > self.pPDiodes
                self.pPDiodes = self.pPDiodes+.1;
                if self.pPDiodes >= 6
                    %Cant increase num diodes anymore
                    self.pPDiodes = 6;
                end
            elseif in == self.pPDiodes
                return;
            else
                self.pPDiodes = self.pPDiodes-.1;
            end
        end
        %-------------------------------------------------------------------
        % 
        %-------------------------------------------------------------------
        function [] = setSeven(self, in)
           if in > self.pNDiodes
                self.pNDiodes = self.pNDiodes+.1;
                if self.pNDiodes >= 6
                    %Cant increase num diodes anymore
                    self.pNDiodes = 6;
                end
            elseif in == self.pNDiodes
                return;
            else
                self.pNDiodes = self.pNDiodes-.1;
            end
        end
        %-------------------------------------------------------------------
        % Set Diodes, helper function, calls each of the harms function
        %-------------------------------------------------------------------
        function [] = setDiodes(self, numHarm, in)
            self.isSetup = false;
            switch numHarm
                case 2
                    setTwo(self, in);
                case 3
                    setThree(self, in);
                case 4
                    setFour(self, in);
                case 5
                    setFive(self, in);
                case 6
                    setSix(self, in);
                case 7
                    setSeven(self, in);
                otherwise
                    return;
            end
        end
        %-------------------------------------------------------------------
        % Calculate local fft peaks
        %-------------------------------------------------------------------
        function [] = calcHarms(self, in)
            exist = false;
                for i = 1:3
                    if self.B(i)
                        exist = true;
                    end
                end
                for i = 1:3
                    if self.A(i)
                        exist = true;
                    end
                end
                    [N, channels] = size(in);
                    Fs = self.getSampleRate();
                    self.sampsPerFrame = N;
                    freq = 0:Fs/N:Fs/2;
                    self.pBins = freq;
                    xdft = fft(in(:, 1));
                    xdft = xdft(1:floor(N/2+1));
                    psdx = (1/(Fs*N)) * abs(xdft).^2;
                    psdx(2:end-1) = 2*psdx(2:end-1);
                    psdx = 10 * log(psdx);
                    %Find peaks in fft
                    if N > 2
                        [x, loc] = findpeaks(psdx, freq);
                        for i = 1:length(loc)
                            self.pHarms(i) = loc(i);
                            if self.pHarms(i) == 0
                               self.pHarms(i) = self.pHarms(i-1) * 2;
                            end
                        end
                        if self.pHarms(self.attenFreq) >= 20e3
                            self.buildFiltCoeffs(self.pHarms(self.attenFreq-2));
                            
                        else
                            self.buildFiltCoeffs(self.pHarms(self.attenFreq));
                        end
                        self.isSetup = true;
                    else
                        self.isSetup = false;
                        self.builtFilter = false;
                    end
                    
                    %plot(freq, psdx)
                    %findpeaks(self.localMax(:, 2), self.localMax(:, 1))
                    
                    
               
            end
    end
end