function [] = progress_in_console(i);

if i == 1 
    disp('   '); 
end
if (i < 10) 
    fprintf('\b\b'); 
elseif (i  < 100) 
    fprintf('\b\b\b'); 
else 
    fprintf('\b\b\b\b');  
end
fprintf('%d %s', i, '');    




% if (i < 10) 
%     if i ==1
%         fprintf(' '); 
%     else
%         fprintf('\b\b\b'); 
%     end
% elseif (i  < 100) 
%     fprintf('\b\b\b\b'); 
% else 
%     fprintf('\b\b\b\b\b');  
% end
% fprintf('%d %s', i, '');
% fprintf('\n');