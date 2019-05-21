classdef inner_class_2
   properties
      p_3 = [];
      p_4 = [];
      parent = [];
   end
      
   methods (Access = ?outer_class)
       function obj = inner_class_2(parent)
            obj.parent = parent;
       end
   end
   
   methods
      function obj = set.p_3(obj,in)
          display('in setter of p_3')  
          obj.p_3 = in;
          obj.parent.flat.p_3 = in;
      end
      function out = get.p_3(obj)
         out = obj.p_3;
      end
      function obj = set.p_4(obj,in)
          display('in setter of p_4')  
          obj.p_4 = in;
          obj.parent.flat.p_4 = in;
      end
      function out = get.p_4(obj)
         out = obj.p_4;
      end
   end
end
