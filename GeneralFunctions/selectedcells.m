function newdata = selectedcells(olddata,selection)
% Take an exisiting multiple cell imaging dataset and keep only the
% selected cells.
newdata=olddata;
newdata.cells=[];
newdata.cells=olddata.cells(selection);
end