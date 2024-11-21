% Extract time and data from the timeseries object
timeData = DI.Time;    % Extract the time vector
valueData = DI.Data;   % Extract the data values

% Create a table for better organization
dataTable = table(timeData, valueData, 'VariableNames', {'Time', 'Data'});

% Save the table to an Excel file
writetable(dataTable, 'DI2.xlsx');