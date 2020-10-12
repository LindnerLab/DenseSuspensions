function [Pparticles] = ParticleSelection(Img, Img_norm, r, SelectionCriteria,verbose)
        Ncriteria = size(SelectionCriteria,2);                              % Determine the number of criteria given
        Centroids = table2array(regionprops('table',Img,Img_norm,'WeightedCentroid'));       % Determine the centroids of all potential particles found
        Nregions = size(Centroids,1);                                       % Determine the number of potential particles found

        if Ncriteria == 0
            if verbose
                fprintf('No selection criteria were given, so all centroids are returned.\n');
            end
            Pparticles = Centroids;
            return
        elseif Nregions == 0
            if verbose
                fprintf('No particles were found.\n');
            end
            Pparticles = [];
            return
        end

        Selection = zeros(Ncriteria,Nregions);
        for i = 1:Ncriteria
            Property = SelectionCriteria(i).Property;                          % Region property to check
            Value = SelectionCriteria(i).Value;                                % Value to check the property against
            Criteria = SelectionCriteria(i).Criteria;                          % Criteria to be a valid particle (either 'greater' or 'smaller')

            P_properties = table2array(regionprops('table',Img,Property));  % Determine the values of the regionproperties to compare with Value

            % Adds a vector with logical values to Selection, where 1 means the
            % regionproperty of that region complies with the given criteria
            % and 0 means that the regionproperty doesnt comply with the
            % criteria.
            if strcmpi(Criteria, 'greater')
                Selection(i,:) = logical([P_properties > Value]);
            elseif strcmpi(Criteria, 'smaller')
                Selection(i,:) = logical([P_properties < Value]);
            end
        end

        % prod(Selection,1) will multiply all logical values for each region,
        % meaning that if a region passes all tests (all logical values are 1),
        % the result will be 1, and thus the Centroid of that region will be
        % added to the list. If a region fails one or more tests, the result
        % will be 0, and thus that region will not be selected as a valid
        % particle.
        % Also, r is subtracted from the x and y position of the peak found.
        % This shift is caused by the position of the 'object' on the mask. The
        % object is centered aro    und (r,r), to limit edge effects, but this
        % shift has to be subtracted, which is done here.
        Pparticles = Centroids(logical(prod(Selection,1)),:) - r;

        if isempty(Pparticles) && verbose
            fprintf('No particles met the requirements');
        end
end