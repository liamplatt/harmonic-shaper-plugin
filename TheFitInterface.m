classdef TheFitInterface < audioPlugin
    properties
        GainHarm1 = 0
        GainHarm2 = 0
        GainHarm3 = 0
        GainHarm4 = 0
        GainHarm5 = 0
        GainHarm6 = 0
        Distortion = 1;
        Input = 1;
        Output = 1;
        Enable = true;
    end
    properties (Constant)
        PluginInterface = audioPluginInterface( ...
            audioPluginParameter('GainHarm1', ...
                'Label','dB', ...
                'Mapping',{'lin',-20,0}, ...
                'Style','vslider', ...                                     %<--
                'Layout',[2,3;4,3], ...                                    %<--
                'DisplayName','Low','DisplayNameLocation','Above'), ...    %<--
            audioPluginParameter('GainHarm2', ...
                'Label','dB', ...
                'Mapping',{'lin',-20,0}, ...
                'Style','vslider', ...                                     %<--
                'Layout',[2,4;4,4], ...                                    %<--
                'DisplayName','Low','DisplayNameLocation','Above'), ...    %<--
            audioPluginParameter('GainHarm3', ...
                'Label','dB', ...
                'Mapping',{'lin',-20,0}, ...
                'Style','vslider', ...                                     %<--
                'Layout',[2,5;4,5], ...                                    %<--
                'DisplayName','Low','DisplayNameLocation','Above'), ...    %<--
            audioPluginParameter('GainHarm4', ...
                'Label','dB', ...
                'Mapping',{'lin',-20,0}, ...
                'Style','vslider', ...                                     %<--
                'Layout',[2,6;4,6], ...                                    %<--
                'DisplayName','Low','DisplayNameLocation','Above'), ...    %<--
            audioPluginParameter('GainHarm5', ...
                'Label','dB', ...
                'Mapping',{'lin',-20,0}, ...
                'Style','vslider', ...                                     %<--
                'Layout',[2,7;4,7], ...                                    %<--
                'DisplayName','Low','DisplayNameLocation','Above'), ...    %<--
            audioPluginParameter('GainHarm6', ...
                'Label','dB', ...
                'Mapping',{'lin',-20,0}, ...
                'Style','vslider', ...                                     %<--
                'Layout',[2,8;4,8], ...                                    %<--
                'DisplayName','Low','DisplayNameLocation','Above'), ...    %<--
            audioPluginParameter('Output', ...
                'Mapping',{'lin',0,2}, ...
                'Style','rotaryknob', ...                                  %<--
                'Layout',[3,10], ...                                        %<--
                'DisplayNameLocation','Above'), ...                        %<--
            audioPluginParameter('Input', ...
                'Mapping',{'lin',0,2}, ...
                'Style','rotaryknob', ...                                  %<--
                'Layout',[3,1], ...                                        %<--
                'DisplayNameLocation','Above'), ...                        %<--
            audioPluginParameter('Distortion', ...
                'Mapping',{'lin',0,2}, ...
                'Style','rotaryknob', ...                                  %<--
                'Layout',[6,5;6,6], ...                                        %<--
                'DisplayNameLocation','Above'), ...                        %<--
            audioPluginParameter('Enable', ...
                'Style','vtoggle', ...                                     %<--
                'Layout',[6,1], ...                                        %<--
                'DisplayNameLocation','None'), ...                         %<--
                ...
            audioPluginGridLayout( ...                                     %<--
                'RowHeight',[20,20,160,20,20,150], ...                        %<--
                'ColumnWidth',[150,20,100,100,100,100,100,100,20,150], ...                    %<--
                'Padding',[10,10,10,10]))                                  %<--
    end
end