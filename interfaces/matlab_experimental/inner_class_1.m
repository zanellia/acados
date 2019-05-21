classdef inner_class_1
   properties
      p_1 = [];
      p_2 = [];
      parent = [];
   end
   
   methods (Access = ?outer_class)
       function obj = inner_class_1(parent)
            obj.parent = parent;
       end
   end
   
   methods 
      function obj = set.p_1(obj,in)
          display('in setter of p_1')  
          obj.p_1 = in;
          obj.parent.flat.p_1 = in;
      end
      function out = get.p_1(obj)
         out = obj.p_1;
      end
      function obj = set.p_2(obj,in)
          display('in setter of p_2')  
          obj.p_2 = in;
          obj.parent.flat.p_2 = in;
      end
      function out = get.p_2(obj)
         out = obj.p_2;
      end
   end
end
