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
        public double Phi;
        public Point(int x, int y, double val){
            this.x = x;
            this.y = y;
            value = val;
        }
    }
    Image img;
    public Image res;
    ArrayList<Point> Points;

    public InterestPoint(Image img){
        Points = new ArrayList<>();
        this.img = img;
    }

    public BufferedImage Moravek(double threshold, int radius, int pointsCount) {
        int inc = radius + 1;
        double[] augmented_img = img.ImageSupplement(img.matrix, inc*2+1, inc*2+1);
        int aug_width = img.getWidth() + inc*2;
        int aug_height = img.getHeight() + inc*2;
        double[] matrix = new double[img.getHeight()*img.getWidth()];

        for (int y = 0; y < img.getHeight(); y++) {
            for (int x = 0; x < img.getWidth(); x++) {
                // 8 направлений
                ArrayList<Double> local_S = new ArrayList<Double>();
                for(int i=0; i <8;i++)
                    local_S.add((double)0);
                for (int u = -radius; u < radius; u++) {
                    for (int v = -radius; v < radius; v++) {
                        double[] directDiff = new double[8];
                        //int pixel = img.getPixel(x + u, y + v);
                        double pixel = augmented_img[(y + inc + v) * aug_width + (x + inc + u)];
                        directDiff[0] = pixel - augmented_img[(y + inc + v - 1) * aug_width + (x + inc + u)];
                        directDiff[1] = pixel - augmented_img[(y + inc + v + 1) * aug_width + (x + inc + u)];
                        directDiff[2] = pixel - augmented_img[(y + inc + v) * aug_width + (x + inc + u + 1)];
                        directDiff[3] = pixel - augmented_img[(y + inc + v - 1) * aug_width + (x + inc + u + 1)];
                        directDiff[4] = pixel - augmented_img[(y + inc + v + 1) * aug_width + (x + inc + u + 1)];
                        directDiff[5] = pixel - augmented_img[(y + inc + v) * aug_width + (x + inc + u - 1)];
                        directDiff[6] = pixel - augmented_img[(y + inc + v - 1) * aug_width + (x + inc + u - 1)];
                        directDiff[7] = pixel - augmented_img[(y + inc + v + 1) * aug_width + (x + inc + u - 1)];

                        for (int i = 0; i < 8; i++) {
                            double value = local_S.get(i) + directDiff[i] * directDiff[i];
                            local_S.set(i, value);
                        }
                    }
                }
                matrix[y*img.getWidth()+x] = Collections.min(local_S);
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

        ArrayList<Point> p_new = new ArrayList<Point>();
        for(int i=0; i<Points.size(); i++){
            if(Points.get(i).value >= threshold){
                p_new.add(Points.get(i));
            }
        }
        Points = p_new;
        Image img_S = new Image(matrix, img.getHeight(), img.getWidth());
        Points = localMaximum(Points, img_S);
        Points = anmsFilter(Points, pointsCount);

        BufferedImage newImage = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);

        Graphics2D g = newImage.createGraphics();
        g.drawImage(img.bufimg, 0, 0, img.getWidth(), img.getHeight(), null);
        BasicStroke bs =new BasicStroke(1);
        g.setStroke(bs);
        int pc=0;
        if (pointsCount > Points.size()) pc = Points.size();
        else pc = pointsCount;
        for(int i=0;i<pc;i++){
            Point p = Points.get(i);
            Ellipse2D.Float circle = new Ellipse2D.Float(p.x, p.y, 3, 3);
            g.setPaint(Color.white);
            g.fill(circle);
            g.setPaint(Color.black);
            g.draw(circle);
        }
        g.dispose();
        return newImage;
    }

    public ArrayList<Point> anmsFilter(ArrayList<Point> points, int pointsCount) {
        Boolean[] flagUsedPoints = new Boolean[points.size()];

        for (int i = 0; i < points.size(); i++)
            flagUsedPoints[i] = true;

        int radius = 0;
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
                for (int j = 0; j < points.size(); j++)
                {
                    if (flagUsedPoints[j])
                    {
                        Point p2 = points.get(j);
                        if ((p1.value > p2.value) && (Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)) <= radius))
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

    public BufferedImage HarrisMap (int radius){
        double[] gauss = img.GaussKernel((double) radius/3);

        double[] image_dx = img.DerivativeX(false);
        double[] image_dy = img.DerivativeY(false);

        int inc = radius + 1;

        double[] aug_img_dx = img.ImageSupplement(image_dx, inc*2+1, inc*2+1);
        double[] aug_img_dy = img.ImageSupplement(image_dy, inc*2+1, inc*2+1);
        double[] matrix = new double[img.getHeight()*img.getWidth()];

        for (int x = 0; x < img.getWidth(); x++)
        {
            for (int y = 0; y < img.getHeight(); y++)
            {
                double val = lambda(aug_img_dx, aug_img_dy, x, y, radius, gauss);
                matrix[y*img.getWidth()+x] = val;
            }
        }

        Image img_S = new Image(matrix, img.getHeight(), img.getWidth());
        img_S.matrix = img_S.NormalizeTo255(img_S.matrix);
        BufferedImage bufimg = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        bufimg.getRaster().setSamples(0,0, img.getWidth(), img.getHeight(),0, img_S.matrix);

        return bufimg;
    }

    public BufferedImage Harris (double threshold, int radius, int pointsCount){
        double[] gauss = img.GaussKernel((double) radius/3);

        double[] image_dx = img.DerivativeX(false);
        double[] image_dy = img.DerivativeY(false);

        int inc = radius + 1;

        double[] aug_img_dx = img.ImageSupplement(image_dx, inc*2+1, inc*2+1);
        double[] aug_img_dy = img.ImageSupplement(image_dy, inc*2+1, inc*2+1);
        double[] matrix = new double[img.getHeight()*img.getWidth()];

        for (int x = 0; x < img.getWidth(); x++)
        {
            for (int y = 0; y < img.getHeight(); y++)
            {
                double val = lambda(aug_img_dx, aug_img_dy, x, y, radius, gauss);
                Points.add(new Point(x, y, val));
                matrix[y*img.getWidth()+x] = val;
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

        ArrayList<Point> p_new = new ArrayList<Point>();
        for(int i=0; i<Points.size(); i++){
            if(Points.get(i).value >= threshold){
                p_new.add(Points.get(i));
            }
        }
        Points = p_new;
        Image img_S = new Image(matrix, img.getHeight(), img.getWidth());

        ArrayList<Point> localMaximumPoints = localMaximum(Points, img_S);
        Points = anmsFilter(localMaximumPoints, pointsCount);

        BufferedImage newImage = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = newImage.createGraphics();
        g.drawImage(img.bufimg, 0, 0, img.getWidth(), img.getHeight(), null);
        BasicStroke bs =new BasicStroke(1);
        g.setStroke(bs);
        int pc=0;
        if (pointsCount > Points.size()) pc = Points.size();
        else pc = pointsCount;
        for(int i=0;i<pc;i++){
            Point p = Points.get(i);
            Ellipse2D.Float circle = new Ellipse2D.Float(p.x-1, p.y-1, 3, 3);
            g.setPaint(Color.white);
            g.fill(circle);
            g.setPaint(Color.black);
            g.draw(circle);
        }
        g.dispose();
        res = new Image(newImage);
        return newImage;
    }

    public double lambda(double[] image_dx, double[] image_dy, int x, int y, int radius, double[] gauss)
    {
        int inc = radius + 1;

        int aug_width = img.getWidth() + inc*2;
        int aug_height = img.getHeight() + inc*2;

        double A = 0, B = 0, C = 0;
        int k = 0, q = 0;
        int g_width = (int)Math.sqrt(gauss.length);
        for (int i = x - radius; i < x + radius; i++)
        {
            for (int j = y - radius; j < y + radius; j++)
            {
                //double curA = image_dx[j*img.getHeight() + i];
                //double curB = image_dy[j*img.getHeight() + i];
                double curA = image_dx[(j+inc)*aug_width + i+inc];
                double curB = image_dy[(j+inc)*aug_width + i+inc];
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

    public ArrayList<Point> localMaximum(ArrayList<Point> points, Image image_S)
    {
        ArrayList<Point> result = new ArrayList<Point>();
        int radius = 2;

        //Смотрим интенсивность всех точек в окрестности
        //c заданным радиусом
        for (int i = 0; i < points.size(); i++)
        {
            Point p1 = points.get(i);
            boolean flagMaximum = true;

            for (int j = -radius; j <= radius; ++j)
            {
                for (int k = -radius; k <= radius; ++k)
                {
                    if (j == 0 && k == 0)
                    {
                        continue;
                    }

                    if (image_S.GetPixel(p1.x + j, p1.y + k) >= p1.value)
                    {
                        flagMaximum = false;
                        break;
                    }
                }
            }

            if (flagMaximum == true)
            {
                result.add(p1);
            }
        }
        return result;
    }
}
