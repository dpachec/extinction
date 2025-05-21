function d = cohens_d(x1, x2)
    % Compute the means of both groups
    M1 = mean(x1);
    M2 = mean(x2);

    % Compute the standard deviations of both groups
    S1 = std(x1);
    S2 = std(x2);

    % Compute the sample sizes
    n1 = length(x1);
    n2 = length(x2);

    % Compute the pooled standard deviation
    Sp = sqrt(((n1 - 1) * S1^2 + (n2 - 1) * S2^2) / (n1 + n2 - 2));

    % Compute Cohen's d
    d = (M1 - M2) / Sp;
end
