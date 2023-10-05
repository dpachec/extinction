function r = nancorr_hui( X,Y,corrtype,crit )
% NANCORR calculates the sample correlation coefficient
%    for the series with NaNs expected.
%    X is the one series, Y is another.
%    chanCrit is the least number of channels left remove NaN
% by Hui Jan.23,2019
% X = test1;
% Y = test2;
X=X(:);
Y=Y(:);
L1=length(X);
L2=length(Y);
if L1 ~= L2
    error('The samples must be of the same length')
end
for i=1:L1
    if isnan(X(i))
        Y(i)=nan;
    end
    if isnan(Y(i))
        X(i)=nan;
    end
end

if sum(~isnan(X)) < crit
    r = NaN;
else
    if strcmp(corrtype,'pearson')
        Xm=nanmean(X);
        Ym=nanmean(Y);
        r=nansum((X-Xm).*(Y-Ym))/sqrt((nansum((X-Xm).^2))*(nansum((Y-Ym).^2)));
    elseif strcmp(corrtype,'spearman')
        X1 = X(~isnan(X));
        Y1 = Y(~isnan(Y));

        % Find the data length
        N = length(X1);
        % Get the ranks of x
        [~,R1] = sort(X1,'descend');
        [~,z2] = sort(R1);
        R2 = (1:length(X1))';
        R_X=R2(z2);

        % Get the ranks of y
        [~,R1] = sort(Y1,'descend');
        [~,z2] = sort(R1);
        R2 = (1:length(Y1))';
        R_Y=R2(z2);

        % Calculate the correlation coefficient
        r = 1-6*sum((R_X-R_Y).^2)/N/(N^2-1);
    end
end

end