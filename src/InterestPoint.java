import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class InterestPoint {

    public  class Point{
        public int x;
        public int y;
        public double value;
        public Point(int x, int y, double val){
            this.x = x;
            this.y = y;
            value = val;
        }
    }
    Image img;
    ArrayList<Point> Points;

    public InterestPoint(Image img){
        Points = new ArrayList<>();
        this.img = img;
    }
    public BufferedImage Moravek(double threshold, int radius, int pointsCount) {
        int inc = radius + 1;
        int[] augmented_img = img.ImageSupplement(img.matrix, inc*2+1, inc*2+1);
        int aug_width = img.getWidth() + inc*2;
        int aug_height = img.getHeight() + inc*2;

        for (int y = 0; y < img.getHeight(); y++) {
            for (int x = 0; x < img.getWidth(); x++) {
                //double[] local_S = new double[8];                  // 8 направлений
                ArrayList<Double> local_S = new ArrayList<Double>();
                for(int i=0; i <8;i++)
                    local_S.add((double)0);
                for (int u = -radius; u < radius; u++) {
                    for (int v = -radius; v < radius; v++) {
                        double[] directDiff = new double[8];
                        //int pixel = img.getPixel(x + u, y + v);
                        int pixel = augmented_img[(y + inc + v) * aug_width + (x + inc + u)];
                        directDiff[0] = (double) (pixel - augmented_img[(y + inc + v - 1) * aug_width + (x + inc + u)])/255;
                        directDiff[1] = (double) (pixel - augmented_img[(y + inc + v + 1) * aug_width + (x + inc + u)])/255;
                        directDiff[2] = (double) (pixel - augmented_img[(y + inc + v) * aug_width + (x + inc + u + 1)])/255;
                        directDiff[3] = (double) (pixel - augmented_img[(y + inc + v - 1) * aug_width + (x + inc + u + 1)])/255;
                        directDiff[4] = (double) (pixel - augmented_img[(y + inc + v + 1) * aug_width + (x + inc + u + 1)])/255;
                        directDiff[5] = (double) (pixel - augmented_img[(y + inc + v) * aug_width + (x + inc + u - 1)])/255;
                        directDiff[6] = (double) (pixel - augmented_img[(y + inc + v - 1) * aug_width + (x + inc + u - 1)])/255;
                        directDiff[7] = (double) (pixel - augmented_img[(y + inc + v + 1) * aug_width + (x + inc + u - 1)])/255;

                        for (int i = 0; i < 8; i++) {
                            double value = local_S.get(i) + directDiff[i] * directDiff[i];
                            local_S.set(i, value);
                        }
                    }
                }
                Points.add(new Point(x, y, Collections.min(local_S)));
            }
        }
        Points.sort(new Comparator<Point>() {
            @Override
            public int compare(Point o1, Point o2) {
                if (o1.value == o2.value) return 0;
                else if (o1.value< o2.value) return 1;
                else return -1;
            }
        });

        for(int i=0; i<Points.size(); i++){
            if(Points.get(i).value < threshold){
                Points.remove(i);
                i--;
            }
        }

        Points = anmsFilter(Points, pointsCount);

        BufferedImage newImage = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);

        Graphics2D g = newImage.createGraphics();
        g.drawImage(img.bufimg, 0, 0, img.getWidth(), img.getHeight(), null);
        float[] shtrich = {14, 5};
        BasicStroke bs =new BasicStroke(2);
        g.setStroke(bs);
        g.setPaint(Color.red);

        for(int i=0;i<Points.size();i++){
            Point p = Points.get(i);
            g.draw(new Ellipse2D.Float(p.x, p.y, 1, 1));
        }
        g.dispose();
        return newImage;
    }

    public ArrayList<Point> anmsFilter(ArrayList<Point> points, int pointsCount) {
        Boolean[] flagUsedPoints = new Boolean[points.size()];

        for (int i = 0; i < points.size(); i++)
            flagUsedPoints[i] = true;

        int radius = 3;
        int usedPointsCount = points.size();
        while (usedPointsCount > pointsCount)
        {
            for (int i = 0; i < points.size(); i++)
            {
                if (!flagUsedPoints[i])
                {
                    continue;
                }

                Point p1 = points.get(i);
                for (int j = i + 1; j < points.size(); j++)
                {
                    if (flagUsedPoints[j])
                    {
                        Point p2 = points.get(j);
                        if (p1.value * 0.9 > p2.value && Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)) <= radius)
                        {
                            flagUsedPoints[j] = false;
                            usedPointsCount--;
                            if (usedPointsCount <= pointsCount)
                            {
                                break;
                            }
                        }
                    }
                }
            }
            radius++;
        }
        ArrayList<Point> resultPoints = new ArrayList<Point>();
        for (int i = 0; i < points.size(); i++)
        {
            if (flagUsedPoints[i])
            {
                resultPoints.add(points.get(i));
            }
        }
        return resultPoints;
    }

    public BufferedImage Harris (double threshold, int radius, int pointsCount){

        /*
        BufferedImage g = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        g.getRaster().setSamples(0,0, img.getWidth(), img.getHeight(),0, img.matrix);
        Image gauss = new Image(g);
        gauss.GaussianBlur(radius/3);
        //Kernel gauss = KernelCreator.getGaussSlowPoke((double)radius / 3);
        //image = ImageConverter.convolution(image, gauss);
        */

        double[] gauss = img.GaussKernel(radius/3);

        int[] image_dx = img.DerivativeX();
        int[] image_dy = img.DerivativeY();
        for (int x = 0; x < img.getWidth(); x++)
        {
            for (int y = 0; y < img.getHeight(); y++)
            {
                Points.add(new Point(x,y, lambda(image_dx, image_dy, x, y, radius, gauss)));
                //image_S.setPixel(x, y, lambda(image_dx, image_dy, x, y, radius, gauss));
            }
        }

        for(int i=0; i<Points.size(); i++){
            if(Points.get(i).value < threshold){
                Points.remove(i);
                i--;
            }
        }

        List<Point> localMaximumPoints = localMaximum(thresholdFilter(image_S, threshold), image_S);
        return anmsFilter(localMaximumPoints, pointsCount);
    }

    public double lambda(int[] image_dx, int[] image_dy, int x, int y, int radius, double[] gauss)
    {
        double A = 0, B = 0, C = 0;
        int k = 0, q = 0;
        int g_width = (int)Math.sqrt(gauss.length);
        for (int i = x - radius; i < x + radius; i++)
        {
            for (int j = y - radius; j < y + radius; j++)
            {
                int curA = image_dx[j*img.getHeight() + i];
                int curB = image_dy[j*img.getHeight() + i];
                A += curA * curA * gauss[q*g_width + k];
                B += curA * curB * gauss[q*g_width + k];
                C += curB * curB * gauss[q*g_width + k];
                k++;
            }
            k = 0;
            q++;
        }
        return ((A * C - B * B) - 0.05 * (A + C) * (A + C)); //вариант оригинального Харриса
    }
}
