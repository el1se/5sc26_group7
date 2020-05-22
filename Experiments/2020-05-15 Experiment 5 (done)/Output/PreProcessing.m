%% changing names
names = Data6.Properties.VariableNames;
names{1,1} = 'x1_mm';
names{1,2} = 'x3_mm';
names{1,3} = 'x2_mm';
names{1,4} = 'x1LeakPos';
names{1,5} = 'x1x3ValvePos';
names{1,6} = 'x2LeakPos';
names{1,7} = 'x3x2ValvePos';
names{1,8} = 'x3LeakPos';
Data6.Properties.VariableNames = names;
clear names