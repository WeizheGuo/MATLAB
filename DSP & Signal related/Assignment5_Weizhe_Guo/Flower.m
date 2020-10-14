classdef Flower
   properties
      petalWidth
      petalLength
      sepalWidth
      sepalLength
      species
   end
    methods
      function flower = Flower(a,b,c,d,A)
         flower.petalWidth = a;
         flower.petalLength = b;
         flower.sepalWidth = c;
         flower.sepalLength = d;
         flower.species = A;   
      end
      function l = getSLength(flower)
         l = flower.sepalLength;
      end
      function report(flower)
          sprintf('The length and width of its sepal are %.1f cm and %.1f cm respectively, while that of its petal are %.1f cm and %.1f cm respectively. It belongs to the %s class. \n'...
              ,flower.petalWidth,flower.petalLength,flower.sepalWidth,flower.sepalLength,flower.species)
      end
   end
end
