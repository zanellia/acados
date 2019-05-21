classdef outer_class < handle
   properties
      a = [];
      b = [];
      flat = [];
   end
   
   methods
       function obj = outer_class()
            obj.a = inner_class_1(obj);
            obj.b = inner_class_2(obj);
       end
   end
end
