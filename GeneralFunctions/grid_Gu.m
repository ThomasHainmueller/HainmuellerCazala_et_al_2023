
n_bin = 80; % has to be 80, equivalent to 5cm/s
n_categories = length(data.metadata.categories);
grid_cell = zeros(n_categories,length(data.cells));

for cell = 1:length(data.cells)    
    for cat = 1:n_categories 
        
        n_trial = length(data.metadata.categories{1,cat}.filename);
        n_total_sample = n_trial*n_bin;
        
        %% data assigning 
        raw_dfoy = zeros(n_trial,n_bin);
        raw_dfoy_concatenated = zeros(1,n_total_sample);

        for k = 1:n_trial
            raw_dfoy(k,:) = data.cells{1,cell}.categories{1,cat}.dFoY{1,k}';
            raw_dfoy_concatenated(n_bin*(k-1)+1:n_bin*k) = raw_dfoy(k,:);
        end    

        raw_dfoy_concatenated(isnan(raw_dfoy_concatenated)) = 0;
        mean_dfoy = mean(raw_dfoy);

        %% dislocated shuffle
        mean_dfoy_shuffled = zeros(1000,n_bin);

        for s = 1:1000

            random_dislocation = randperm(n_total_sample,1);

            while random_dislocation > 0.95*n_total_sample || random_dislocation < 0.05*n_total_sample
                  random_dislocation = randperm(n_total_sample,1);
            end      

            raw_dfoy_concatenated_dislocated = zeros(1,n_total_sample);

            for i = 1:random_dislocation
                raw_dfoy_concatenated_dislocated(i) = raw_dfoy_concatenated(n_total_sample - random_dislocation+i); 
            end          
            for i = 1:(n_total_sample - random_dislocation) 
                raw_dfoy_concatenated_dislocated(i+random_dislocation) = raw_dfoy_concatenated(i);
            end

            dislocated_dfoy = zeros(length(data.cells{1,cell}.categories{1,cat}.dFoY), n_bin);

            for kkk = 1:length(data.cells{1,cell}.categories{1,cat}.dFoY)
                dislocated_dfoy(kkk,:) = raw_dfoy_concatenated_dislocated(((kkk-1)*n_bin+1):kkk*n_bin); 
            end

            mean_dfoy_shuffled(s,:) = mean(dislocated_dfoy) ;

        end  

        p_value = zeros(1,n_bin);

        for t = 1:n_bin
            p_value(t) = length(find(mean_dfoy_shuffled(:,t)> mean_dfoy(t)))/1000;       
        end  

        %% in and out-of field
        infield = zeros(1,n_bin);
        outfield = zeros(1,n_bin);

        if (1-p_value(1))>0.8 && (1-p_value(2))>0.8
           if any(raw_dfoy(:,1)) == 1 || any(raw_dfoy(:,2)) == 1 
              infield(1)=1; infield(2)=1;
           end
        end   

        if (1-p_value(n_bin))>0.8 && (1-p_value(n_bin-1))>0.8
           if any(raw_dfoy(:,n_bin)) == 1 || any(raw_dfoy(:,n_bin-1)) == 1 
           infield(n_bin)=1; infield(n_bin-1)=1;
           end
        end  

        for n = 3:(n_bin-2)        
            if (1-p_value(n))>0.8 && (1-p_value(n+1))>0.8 && (1-p_value(n+2))>0.8
               if any(raw_dfoy(:,n)) == 1 || any(raw_dfoy(:,n+1)) == 1 || any(raw_dfoy(:,n+2)) == 1
                  infield(n)=1;infield(n+1)=1;infield(n+2)=1;
               end   
            end
        end

        for n = 1:n_bin        
            if (1-p_value(n))<0.25 && (1-p_value(n+1))<0.25
               outfield(n)=1;outfield(n+1)=1;
            end
        end

        %% 1D grid classifier
        n_place_field = 0;   

        if infield(1) == 1
           n_place_field = 1;
        end

        n_place_field = n_place_field + length(find(diff(infield)==1));

        if any(n_place_field)
           w = length(find(infield==1))/n_place_field;
           transition_threshold = ceil(n_bin/5/w);
        else
           w = [];
        end   

        max_width = max(diff(find(diff(infield))));

        if infield(1)==1
           width_first_field = find(diff(infield)==-1);
           max_width = max([max_width,width_first_field]);
        end   

        if infield(n_bin)==1
           width_last_field = n_bin -find(diff(infield)==1,1,'last');
           max_width = max([max_width,width_last_field]);
        end

        if infield(1)==1 && infield(n_bin)==1
           max_width = max([max_width,width_first_field,width_last_field]);
        end

        n_Transition = length(find(abs(diff(nonzeros(infield-outfield)))==2));
        n_bins_in_out_field = length(find(infield==1 | outfield==1));    
        mean_dfoy_infield = mean(mean_dfoy(infield == 1));
        mean_dfoy_outfield = mean(mean_dfoy(outfield == 1));

        %% classifying
        if n_place_field>1 && n_Transition>transition_threshold && max_width<5*w && ...
           n_bins_in_out_field>0.3*n_bin && mean_dfoy_infield./mean_dfoy_outfield >2        
           grid_cell(cat,cell) =1;
        else
           grid_cell(cat,cell) =0;       
        end               
    end 
end
    