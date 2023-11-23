function result=placefield_evaluator(tdata)
%dFabsThreshold=[0.08,0.1,0.15];
%dFrelThreshold = [2,3,5];
%Ftransient=[0.2,0.25,0.3];
%smooth=[true,false]; Parameter set 160610
dFabsThreshold=[0.0,0.05,0.1,0.15];
dFrelThreshold = [3,5,7,9];
Ftransient=[0.2,0.25,0.3];
smooth=[true];
% 160612 Ideal Params 0.0; 7; 0.2; true based on 160501 dataset

for abs = 1:length(dFabsThreshold)
    for rel = 1:1:length(dFrelThreshold)
        for trans=1:length(Ftransient)
            for aver=1
                fprintf('Processed absTsh %1$f, relTsh %2$f, Ftrans %3$f\n',...
                    dFabsThreshold(abs), dFrelThreshold(rel), Ftransient(trans), smooth(aver));
                resultdata=placefields(tdata,'Dombeck','dFabsThresh',dFabsThreshold(abs),...
                    'dFrelThresh',dFrelThreshold(rel),'Ftransient',Ftransient(trans),...
                    'smooth',smooth(aver),'nshuffle',500,'minrate',1/60);
                result{abs,rel,trans,aver}.Pvalues=resultdata.metadata.Placefield_P;
                for category=1:4
                    result{abs,rel,trans,aver}.PFcells(category,1)=...
                        length(find(resultdata.metadata.Placefield_P(:,category)<0.05));
                    result{abs,rel,trans,aver}.PFcells(category,2)=...
                        length(find(resultdata.metadata.Placefield_P(:,category)<1));
                    x=abs+(length(dFrelThreshold)*(rel-1));
                    y=trans+(length(Ftransient)*(aver-1));
                    result{10,10,10,10}.SignSumar(x,y,category)=...
                        length(find(resultdata.metadata.Placefield_P(:,category)<0.05));
                    result{10,10,10,10}.AllSumar(x,y,category)=...
                        length(find(resultdata.metadata.Placefield_P(:,category)<1));
                end
            end
        end
    end
end
end