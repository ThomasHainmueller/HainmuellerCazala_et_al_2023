start_cell = 1;
end_cell = length(data.cells);

for ii = start_cell:end_cell
        cell_number = ii,
        plot_placeactivityedu(data,ii)
end