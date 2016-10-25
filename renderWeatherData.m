function [ output_args ] = renderWeatherData( input_args )
close all
cvsM = csvread('OUT.csv');

varS={'u' 'v' 'w' 'PI0' 'T0' 'rv0' 'rl0'};
heigS={'Low' 'Mid' 'High'};
for varN=1:size(varS,2)
    figH=figure(varN);,
    set(figH,'Position',[150,100,1800,800]);
    set(figH,'name',[varS{varN}]);
    for height=1:size(heigS,2)
       M= cvsM(height:3:size(cvsM,1),:);%height 1=L 2=M 3=H
        subplot(3,2,(height-1)*2+1);%left
        val=M(:,2+3*(varN-1));
        plot(val)
        title([varS{varN} ' ' heigS{height}]);
        
        subplot(3,2,(height-1)*2+2);%right
        
         val1=M(:,3+3*(varN-1));
         val2=M(:,4+3*(varN-1));
         
         hold on;
         plot(val1);
         plot(val2,'r');
         legend('gridR','Prime');
         hold off;


    end
    saveas(varN,[varS{varN}], 'png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% render T
figH=figure(size(varS,2)+1);
set(figH,'Position',[150,100,1800,800]);
intValues=csvread('initState.csv');
gridZ=52;
hSamples=[2 gridZ/2+1 gridZ-2]
for height=1:size(heigS,2)
    M= cvsM(height:3:size(cvsM,1),:);%height 1=L 2=M 3=H
    subplot(3,1,(height-1)+1);%left
    varN=4;%PI0
    valPI=M(:,2+3*(varN-1));
    varN=5;%T0
    valT0=M(:,2+3*(varN-1));
    val=(valPI+intValues(hSamples(height),7)).*valT0;
    
    plot(val)
    title(['T ' heigS{height}]);
end
saveas(size(varS,2)+1,'T', 'png')

end

