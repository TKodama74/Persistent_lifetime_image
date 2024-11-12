%%Cal_PLTImage.m
%{
    <Input variables>
        img : A 2D matrix of a preprocessed 2D image 
        

%}
function [PLTImg, PLTMap, CentPlot] = Cal_PLTImage(img, roi, qbit, k, sigma)
    %{
    <Input variables>
        img : 2D matrix of a preprocessed 2D image 
        roi : Mask image of a region of interest on the img
        qbit : bit depth of input image
        k : Dimension of topological components (k = 0 : connected component; k = 1 : hole component)
        sigma : Sigma of Gaussian filter (Weight of scalarization)
    <Output variables>
        PLTImg : Persistent lifetime image array for k-th topological component
        PLTMap : Persistent lifetime map array for k-th topological component
        CentPlot : Centroid plot array for k-th topological component
    %}
    disp(['Start calculation for ',num2str(k),'-th topological component']);
    %%Step 1: Calculation of bianary images of a filtration
    img = img.*roi;
    [r_img, c_img] = size(img);
    binary_imgs = zeros([r_img, c_img]);
    for t = 0 : 2^qbit-1                        %Threshold value t
        binary_imgs(:,:,t+1) = img >= t;        %Binary images with pixel values of 1 (>= t) and 0 (< t)
    end

    %%Step 2: Calculation of centroid coodinates and connections of the
    %%k-th persistent components in a filtration of binary iamges
    CompsInfo = persistent_component_analysis(img, binary_imgs, qbit, k);

    %%Step3: Calculation of persistent lifetime images
    [PLTImg, PLTMap, CentPlot] = PLTImage_Calculation(img, CompsInfo, qbit, sigma);
end


function CompsInfo = persistent_component_analysis(img, binary_imgs, qbit, k)
    CompsInfo = cell(1,2^qbit);
    direction = 'decrease'; %Direction of filtration 
    for ii = 1:2^qbit   
        if direction == 'increase'
            binary_img = binary_imgs(:,:,ii); 
        elseif direction == 'decrease'
            binary_img = binary_imgs(:,:,2^qbit-ii+1); 
        end
        binary_img = binary_imgs(:,:,ii);      %Select binary image at threshold t in the filtration
        %Detect k-th components on binary_img
        if k == 0
            connective = 8;
            obj = bwconncomp(binary_img, connective);
            prop = regionprops(obj, "Area","Centroid","PixelList");
        elseif k ==1
            connective = 4;
            invbinary_img = imcomplement(binary_img);
            cc = bwconncomp(invbinary_img, connective);
            labeled = labelmatrix(cc);
            n_BG = labeled(1,1);
            holecc = invbinary_img; % remove surrounding backgrounds
            [x,y] = size(binary_img);
            for xx = 1:x
                for yy = 1:y
                    if labeled(xx,yy) == 0
                        holecc(xx,yy) = 0;
                    elseif labeled(xx,yy) == n_BG
                        holecc(xx,yy) = 0;
                    end
                end
            end
            if sum(holecc,'all') == 0
                continue
            end   
            obj = bwconncomp(holecc, 4);
            prop = regionprops(obj, "Area","Centroid","PixelList");
        end

        %%Calculation of component information at threshold t
        n_comps = length(prop);
        if n_comps == 0
            continue
        end
        CompsInfo{1,ii} = cell(1,n_comps);
        for jj = 1:n_comps
            CompsInfo{1,ii}{1,jj} = cell(1,4); CompsInfo{1,ii}{1,jj}{1,1} = zeros(1,2^qbit); %Connection of component
            CompArea = prop(jj).Area; CompsInfo{1,ii}{1,jj}{1,2} = CompArea; %Area of component
            CompCentroid = prop(jj).Centroid; CompsInfo{1,ii}{1,jj}{1,3} = CompCentroid; %Centroid of component
            CompPixList = prop(jj).PixelList; CompsInfo{1,ii}{1,jj}{1,4} = CompPixList; %pixels of component
        end
        %%Detecting the connection of components
        if isempty(CompsInfo{1,ii}) == 1
            continue
        elseif ii == 1 || isempty(CompsInfo{1,ii-1}) == 1
            for ll = 1:length(CompsInfo{1,ii})
                CompsInfo{1,ii}{1,ll}{1,1}(1,ii) =ll;
            end 
        else
            for ll = 1:length(CompsInfo{1,ii}) %ll : current component label
                %Current component pix list
                cCompPixList = CompsInfo{1,ii}{1,ll}{1,4};
                cCompCentPix = CompsInfo{1,ii}{1,ll}{1,3};
                %Previous component info
                %Concat_label
                Conn_label_vec = zeros(1,length(CompsInfo{1,ii-1}));
                for mm = 1:length(CompsInfo{1,ii-1}) %mm : previous component label
                    %Previous component pixels
                    pCompPixList = CompsInfo{1,ii-1}{1,mm}{1,4};
                    %Detective area
                    DetectiveArea = pCompPixList;
                    %Search connection
                    [r_cPixList, ~] = size(cCompPixList);
                    [r_da, ~] = size(DetectiveArea);
                    for pp = 1:r_cPixList
                        if Conn_label_vec(1,mm) == 1
                            continue
                        end
                        for qq = 1:r_da
                            if cCompPixList(pp,:) == DetectiveArea(qq,:)
                                Conn_label_vec(1,mm) = 1;
                            end
                        end
                    end
                end     
                if sum(Conn_label_vec) == 0
                    CompsInfo{1,ii}{1,ll}{1,1}(1,ii) = ll;
                elseif sum(Conn_label_vec) == 1
                    CompsInfo{1,ii}{1,ll}{1,1} = CompsInfo{1,ii-1}{1,find(Conn_label_vec)}{1,1};
                    CompsInfo{1,ii}{1,ll}{1,1}(1,ii) = ll;
                else
                    overlaps = find(Conn_label_vec);
                    p_lifetime_vec = zeros(1,length(Conn_label_vec));
                    p_area_vec = zeros(1,length(Conn_label_vec));
                    p_cent_mat = zeros(length(Conn_label_vec),2);
                    p_minpix_mat = zeros(length(Conn_label_vec),2);
                    for rr = overlaps
                        p_lifetime_vec(1,rr) = nnz(CompsInfo{1,ii-1}{1,rr}{1,1});
                        p_area_vec(1,rr) = CompsInfo{1,ii-1}{1,rr}{1,2};
                        p_cent_mat(rr,:) = CompsInfo{1,ii-1}{1,rr}{1,3};
                        pCompPixList = sortrows(CompsInfo{1,ii-1}{1,rr}{1,4}, [2 1]);
                        p_minpix_mat(rr,:) = pCompPixList(1,:);
                    end
                    distance_cents_vec=zeros(1,length(Conn_label_vec));
                    for ss = 1: length(Conn_label_vec)
                        if sum(p_cent_mat(ss,:))==0
                            distance_cents_vec(ss) = Inf; 
                            p_minpix_mat(ss,:) =[Inf, Inf];
                        else
                            distance_cents_vec(ss) = pdist2(cCompCentPix, p_cent_mat(ss,:)); 
                        end
                    end
                    l_maxlifetime = find(p_lifetime_vec == max(p_lifetime_vec));
                    l_maxarea = find(p_area_vec == max(p_area_vec));
                    l_mincentdist = find(distance_cents_vec == min(distance_cents_vec));
                    [~,idx_mincood]=sortrows(p_minpix_mat, [2 1]);
                    
                    %%Inheritance criteria for connection of components
                    %%1. Longer prtial lifetime
                    %%2. Larger area
                    %%3. Shorter distance between centroids of current and previous components
                    %%4. Smaller coordinate (x, y) of a contained pixel in the component 
                    if length(l_maxlifetime) ==1
                        CompsInfo{1,ii}{1,ll}{1,1} = CompsInfo{1,ii-1}{1,l_maxlifetime}{1,1};
                    elseif length(l_maxarea) ==1
                        CompsInfo{1,ii}{1,ll}{1,1} = CompsInfo{1,ii-1}{1,l_maxarea}{1,1};
                    elseif  length(l_mincentdist) ==1
                        CompsInfo{1,ii}{1,ll}{1,1} = CompsInfo{1,ii-1}{1,l_mincentdist}{1,1};
                    else
                        CompsInfo{1,ii}{1,ll}{1,1} = CompsInfo{1,ii-1}{1,idx_mincood(1)}{1,1};
                    end
                    CompsInfo{1,ii}{1,ll}{1,1}(1,ii) = ll;
                end
            end
        end
    end
  
    %%Sepalate
    for ii = 2:2^qbit
        %%Connection info at a threshold
        if length(CompsInfo{1,ii}) <=1
            continue
        end
        %%If no. of components are larger than 2, determine the connections
        for jj = 1: length(CompsInfo{1,ii})
            if length(CompsInfo{1,ii}) ==jj
                continue
            end
            
            compcode_base = CompsInfo{1,ii}{1,jj}{1,1};
            p_compcode_base = compcode_base;
            p_compcode_base(ii) = 0;
            if sum(p_compcode_base) ==0
                continue
            end
            for kk = jj+1:length(CompsInfo{1,ii})
                compcode_via = CompsInfo{1,ii}{1,kk}{1,1};
                pcompcode_via = compcode_via;
                pcompcode_via(1,ii) = 0;
                if isequal(p_compcode_base, pcompcode_via) == 0
                    continue
                else
                    %%Inheritance criteria of life time on spalation
                    %%Condition1. lower minimum pixel value in the component
                    %%Condition2. Larger area
                    %%Condition3. Shorter distance between centroids of current and previous components
                    %%Condition4. Smaller coordinate (x, y) of a contained pixel in the component 
                    %Condition parameters for base component
                    basecomp_area = CompsInfo{1,ii}{1,jj}{1,2};
                    basecomp_centcood = CompsInfo{1,ii}{1,jj}{1,3};
                    basecomp_PixList = sortrows( CompsInfo{1,ii}{1,jj}{1,4}, [2 1]);
                    basecomp_minpix = basecomp_PixList(1,:);

                    basemask = Inf(size(img));
                    [nr_basecomp_PixList,~] = size(basecomp_PixList);
                    for r_base = 1:nr_basecomp_PixList
                        basecood = basecomp_PixList(r_base,:);
                        basemask(basecood(2),basecood(1)) = 0;
                    end
                    baseimgmasked = double(img) + basemask;
                    basecomp_minval = min(baseimgmasked,[],'all');

                    %Condition parameters for via component
                    viacomp_area = CompsInfo{1,ii}{1,kk}{1,2};
                    viacomp_centcood = CompsInfo{1,ii}{1,kk}{1,3};
                    viacomp_PixList = sortrows(CompsInfo{1,ii}{1,kk}{1,4}, [2 1]);
                    viacomp_minpix = viacomp_PixList(1,:);

                    viamask = Inf(size(img));
                    [nr_viacomp_PixList,~] = size(viacomp_PixList);
                    for r_via = 1:nr_viacomp_PixList
                        viacood = viacomp_PixList(r_via,:);
                        viamask(viacood(2),viacood(1)) = 0;
                    end
                    viaimgmasked = double(img) + viamask;
                    viacomp_minval = min(viaimgmasked,[],'all');

                    for ll  = 1: length(CompsInfo{1,ii-1})
                        if isequal(p_compcode_base, CompsInfo{1,ii-1}{1,ll}{1,1}) == 1
                            precomp_centcood = CompsInfo{1,ii-1}{1,ll}{1,3};
                        end
                    end
                    %Calculate centroid distance from base or via comps to previous comp
                    basecomp_centdist = pdist2(precomp_centcood, basecomp_centcood); 
                    viacomp_centdist = pdist2(viacomp_centcood, basecomp_centcood); 

                    %%Inheritation of lifetime
                    if basecomp_minval~=viacomp_minval
                        if basecomp_minval < viacomp_minval
                            %CompsInfo{1,ii}{1,jj}{1,1}: keep base change via
                            CompsInfo{1,ii}{1,kk}{1,1}(1:ii-1) = 0;
                        else
                            CompsInfo{1,ii}{1,jj}{1,1}(1:ii-1) = 0;
                        end
                    elseif basecomp_area~=viacomp_area
                        if basecomp_area > viacomp_area
                            %CompsInfo{1,ii}{1,jj}{1,1}: keep base change via
                            CompsInfo{1,ii}{1,kk}{1,1}(1:ii-1) = 0;
                        else
                            CompsInfo{1,ii}{1,jj}{1,1}(1:ii-1) = 0;
                        end
                    elseif basecomp_centdist~=viacomp_centdist
                        if basecomp_centdist < viacomp_centdist
                            %CompsInfo{1,ii}{1,jj}{1,1}: keep base change via
                            CompsInfo{1,ii}{1,kk}{1,1}(1:ii-1) = 0;
                        else
                            CompsInfo{1,ii}{1,jj}{1,1}(1:ii-1) = 0;
                        end
                    elseif basecomp_minpix(2)~= viacomp_minpix(2)
                        if basecomp_minpix(2) < viacomp_minpix(2)
                            %CompsInfo{1,ii}{1,jj}{1,1}: keep base change via
                            CompsInfo{1,ii}{1,kk}{1,1}(1:ii-1) = 0;
                        else
                            CompsInfo{1,ii}{1,jj}{1,1}(1:ii-1) = 0;
                        end
                    elseif basecomp_minpix(1)~= viacomp_minpix(1)
                        if basecomp_minpix(1) < viacomp_minpix(1)
                            %CompsInfo{1,ii}{1,jj}{1,1}: keep base change via
                            CompsInfo{1,ii}{1,kk}{1,1}(1:ii-1) = 0;
                        else
                            CompsInfo{1,ii}{1,jj}{1,1}(1:ii-1) = 0;
                        end
                    end
                end
            end
        end
    end
    
    if direction == 'decrease'
        CompsInfo = flip(CompsInfo);
    end
end


function [PLTImg, PLTMap, CentPlot] = PLTImage_Calculation(img, CompsInfo, qbit, sigma)
    imgsize = size(img);
    CentPlot = zeros(imgsize(1),imgsize(2),2^qbit);
    PLTMap = zeros(imgsize(1),imgsize(2),2^qbit);
    PLTImg = zeros(imgsize(1),imgsize(2),2^qbit);

    Msg = ['Calculation of PLT image for sigma_',num2str(sigma)];
    for t = 0:2^qbit-1
        %%Progress bar
        textwaitbar(t, 2^qbit-1, Msg);
        pause(0.005);
        %%Centroid coodinates of persistent components at threshold t
        if isempty(CompsInfo{1,t+1}) ~= 1
          for comp_label = 1:length(CompsInfo{1,t+1})
              cent_cood  = CompsInfo{1,t+1}{1,comp_label}{1,3};
              CentPlot(round(cent_cood(2)),round(cent_cood(1)),t+1) = 1;
          end
        end

        %%Calculation of PLT
        if isempty(CompsInfo{1,t+1}) ~= 1 
          for comp_label = 1:length(CompsInfo{1,t+1})  %%aa: Label of conponemt
              lifetime = nnz(CompsInfo{1,t+1}{1,comp_label}{1,1});
              cent_cood  = round(CompsInfo{1,t+1}{1,comp_label}{1,3});
              PLTMap(cent_cood(2),cent_cood(1),t+1) = lifetime;

              %Gaussian scalarization
              plt_scalar = zeros(imgsize(1),imgsize(2));
              for x = 1:imgsize(1)
                  for y = 1:imgsize(2)
                      plt_scalar(x,y) = lifetime*exp(  -( ((x-cent_cood(2))^2) + ((y-cent_cood(1))^2) ) / (2*(sigma^2)) );
                  end
              end
              PLTImg(:,:,t+1) = PLTImg(:,:,t+1) + plt_scalar;
          end
        end  
    end
end
