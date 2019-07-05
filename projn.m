function z=projn(x)
if x>=0
    z=x;
else
    z=0;
end
end
% function z=projc(x,y,u)
% if y==0
%     z=0;
% else
% if abs(x)<=u*abs(y)
%     z=x/abs(y);
% else
%     if x>u*abs(y)
%     z=u;
%     else 
%         z=-u;
%     end
% end
% end
% end