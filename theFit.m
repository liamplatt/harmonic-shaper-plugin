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
    end
    %-----------------------------------------------------------------------
    % Public Properties - End user interacts with these
    %-----------------------------------------------------------------------
    properties
        PDiodes = 1;
        NDiodes = 1;
        
    end
    %-----------------------------------------------------------------------
    % Constant Properties - Used to define plugin interface
    %-----------------------------------------------------------------------
    properties (Constant)
        PluginInterface = audioPluginInterface( ...
            audioPluginParameter('PDiodes', ...
            'Label','#', ...
            'Mapping',{'lin',0,4}),...
            audioPluginParameter('NDiodes', ...
            'Label','#', ...
            'Mapping',{'lin',0,4}));
    end
    
    methods
        %-------------------------------------------------------------------
        % Main processing function
        %-------------------------------------------------------------------
        function y = process(self,x)
            [N, channels] = size(x)
            
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
        
        function set.Harmonics(self, in)
            
            
        end
        
    end
    methods (Access = private)
        %-------------------------------------------------------------------
        % Calculate Filter Coefficients
        %-------------------------------------------------------------------
        function [B,A] = highpassCoeffs(~,fc,fs)
            
        end
    end
end