function result = HCPIdx (input, whtIndx)
switch whtIndx
    case 'dir'
        len = length(input);
        if (len == 3)
            UU = input(1);
            VV = input(2);
            WW = input(3);
            uu = (2*UU-VV)/3;
            vv = (2*VV-UU)/3;
            tt = -(uu+vv);
            ww = WW;
            result = [uu,vv,tt,ww];
        elseif (len == 4)
            uu = input(1);
            vv = input(2);
            ww = input(4);
            UU = 2*uu + vv;
            VV = uu + 2*vv;
            WW = ww;
            result = [UU,VV,WW];
        else
            disp('dir input array should have a length of 4 or 3');
        end
      case 'face'
         len = length(input);
         if (len == 3)
            UU = input(1);
            VV = input(2);
            WW = input(3);
            uu = UU;
            vv = VV;
            tt = -(UU+VV);
            ww = WW;
            result = [uu,vv,tt,ww];
         elseif (len == 4)
            uu = input(1);
            vv = input(2);
            ww = input(4);
            UU = uu;
            VV = vv;
            WW = ww;
            result = [UU,VV,WW];
         else
            disp('dir input array should have a length of 4 or 3');
         end
end
end