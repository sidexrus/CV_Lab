import java.awt.image.BufferedImage;
import java.lang.management.BufferPoolMXBean;
import java.util.ArrayList;

public class DescriptorCreator {

    public static class Vector
    {
        public Descriptor first;
        public Descriptor second;

        public Vector(Descriptor first, Descriptor second)
        {
            this.first = first;
            this.second = second;
        }
    }

    public static Descriptor[] getDescriptors(Image image, ArrayList<InterestPoint.Point> interestPoints, int radius, int basketCount, int barCharCount) {
        int dimension = 2 * radius; //размерность окрестности интересной точки
        double sector = 2 * Math.PI / basketCount; //размер одной корзины в гистограмме
        double halfSector = Math.PI / basketCount; // размер половины одной корзины в гистограмме
        int barCharStep = dimension / (barCharCount / 4); //шаг гистограммы
        int barCharCountInLine = (barCharCount / 4);


        double[] image_dx = image.DerivativeX(false);
        double[] image_dy = image.DerivativeY(false);
        int aug_width = image.getWidth() + radius*2;
        int aug_height = image.getHeight() + radius*2;
        double[] aug_img_dx = image.ImageSupplement(image_dx, radius*2+1, radius*2+1);
        double[] aug_img_dy = image.ImageSupplement(image_dy, radius*2+1, radius*2+1);

        Descriptor[] descriptors = new Descriptor[interestPoints.size()];
        for (int k = 0; k < interestPoints.size(); k++)
        {
            descriptors[k] = new Descriptor(barCharCount * basketCount, interestPoints.get(k));

            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    int x = interestPoints.get(k).x + radius;
                    int y = interestPoints.get(k).y + radius;
                    double gradient_X = aug_img_dx[(j - radius + y)*aug_width + i - radius + x];
                    double gradient_Y = aug_img_dy[(j - radius + y)*aug_width + i - radius + x];
                    // get value and phi
                    double value = getGradientValue(gradient_X, gradient_Y);
                    double phi = getGradientDirection(gradient_X, gradient_Y);

                    // получаем индекс корзины в которую входит phi и смежную с ней
                    int firstBasketIndex = (int)Math.floor(phi / sector);
                    int secondBasketIndex = (int)(Math.floor((phi - halfSector) / sector) + basketCount) % basketCount;

                    // вычисляем центр
                    double mainBasketPhi = firstBasketIndex * sector + halfSector;

                    // распределяем L(value)
                    double mainBasketValue = (1 - (Math.abs(phi - mainBasketPhi) / sector)) * value;
                    double sideBasketValue = value - mainBasketValue;

                    // вычисляем индекс куда записывать значения
                    int tmp_i = i / barCharStep;
                    int tmp_j = j / barCharStep;

                    int indexMain = (tmp_i * barCharCountInLine + tmp_j) * basketCount + firstBasketIndex;
                    int indexSide = (tmp_i * barCharCountInLine + tmp_j) * basketCount + secondBasketIndex;

                    if (indexMain >= descriptors[k].data.length)
                        indexMain = 0;

                    if (indexSide >= descriptors[k].data.length)
                        indexSide = 0;

                    // записываем значения
                    descriptors[k].data[indexMain] += mainBasketValue;
                    descriptors[k].data[indexSide] += sideBasketValue;
                }
            }
            descriptors[k].normalize();
            descriptors[k].clampData(0, 0.2);
            descriptors[k].normalize();
        }
        return descriptors;
    }

    private static double getGradientValue(double x, double y) { return Math.sqrt(x * x + y * y);    }

    //угол градиента в радианах
    private static double getGradientDirection(double x, double y) { return Math.atan2(y, x) + Math.PI; }

    /* Ориентация точки */
    private static double[] getPointOrientation(double[] image_dx, double[] image_dy, InterestPoint.Point point, int sigma, int radius, int w) {
        int basketCount = 36;

        int dimension = radius * 2;
        double sector = 2 * Math.PI / basketCount;
        double halfSector = Math.PI / basketCount;

        double[] baskets = new double[basketCount];

        for (int i = 0; i < basketCount; i++)
            baskets[i] = 0;

        for (int i = 1; i < dimension; i++)
        {
            for (int j = 1; j < dimension; j++)
            {

                int x = point.x + radius;
                int y = point.y + radius;
                double gradient_X = image_dx[(j - radius + y)*w + i - radius + x];
                double gradient_Y = image_dy[(j - radius + y)*w + i - radius + x];

                // получаем значение(домноженное на Гаусса) и угол
                double value = getGradientValue(gradient_X, gradient_Y)/* * KernelCreator::getGaussValue(i, j, sigma*2, radius)*/;
                double phi = getGradientDirection(gradient_X, gradient_Y);

                // получаем индекс корзины в которую входит phi и смежную с ней
                int firstBasketIndex = (int)Math.floor(phi / sector);
                int secondBasketIndex = (int)(Math.floor((phi - halfSector) / sector) + basketCount) % basketCount;

                // вычисляем центр
                double mainBasketPhi = firstBasketIndex * sector + halfSector;

                // распределяем L(value)
                double mainBasketValue = (1 - (Math.abs(phi - mainBasketPhi) / sector)) * value;
                double sideBasketValue = value - mainBasketValue;

                // записываем значения
                firstBasketIndex = (int) Descriptor.Clamp(0, basketCount - 1, firstBasketIndex);
                secondBasketIndex = (int) Descriptor.Clamp(0, basketCount - 1, secondBasketIndex);
                baskets[firstBasketIndex] += mainBasketValue;
                baskets[secondBasketIndex] += sideBasketValue;
            }
        }

        // Ищем Пики
        double peak_1 = getPeak(baskets, -1);
        double peak_2 = getPeak(baskets, (int)peak_1);

        //хотя бы peak_1 должен быть!

        if (peak_2 != -1 && baskets[(int)peak_2] / baskets[(int)peak_1] >= 0.8)
        { // Если второй пик не ниже 80%
            double[] peaks = new double[2];
            peaks[0] = parabaloidInterpolation(baskets, (int)peak_1);
            peaks[1] = parabaloidInterpolation(baskets, (int)peak_2);
            return peaks;
        }
        else
        {
            double[] peaks = new double[1];
            peaks[0] = parabaloidInterpolation(baskets, (int)peak_1);
            return peaks;
        }
    }

    /* Поиск пика */
    private static double getPeak(double[] baskets, int notEqual) {
        int maxBasketIndex = -1;
        for (int i = 0; i < baskets.length; i++)
        {
            if (baskets[i] > baskets[(i - 1 + baskets.length) % baskets.length]
                    && baskets[i] > baskets[(i + 1) % baskets.length] && i != notEqual)
            {
                if (maxBasketIndex != -1 && baskets[maxBasketIndex] > baskets[i])
                {
                    continue;
                }
                maxBasketIndex = i;
            }
        }
        return maxBasketIndex;
    }

    /* Интерполяция параболой */
    private static double parabaloidInterpolation(double[] baskets, int maxIndex)
    {
        // берём левую и правую корзину и интерполируем параболой
        double left = baskets[(maxIndex - 1 + baskets.length) % baskets.length];
        double right = baskets[(maxIndex + 1) % baskets.length];
        double mid = baskets[maxIndex];

        double sector = 2 * Math.PI / baskets.length;
        double phi = (left - right) / (2 * (left + right - 2 * mid));
        return (phi + maxIndex) * sector + (sector / 2);
    }

    /*  Инвариантость к вращению  */
    public static Descriptor[] getDescriptorsInvRotation(Image image, ArrayList<InterestPoint.Point> interestPoints,
                                                         int radius, int basketCount, int barCharCount) {
        int sigma = 20;
        int dimension = 2 * radius;
        double sector = 2 * Math.PI / basketCount;
        double halfSector = Math.PI / basketCount;
        int barCharStep = dimension / (barCharCount / 4);
        int barCharCountInLine = (barCharCount / 4);

        double[] image_dx = image.DerivativeX(false);
        double[] image_dy = image.DerivativeY(false);
        int aug_width = image.getWidth() + radius*2;
        int aug_height = image.getHeight() + radius*2;
        double[] aug_img_dx = image.ImageSupplement(image_dx, radius*2+1, radius*2+1);
        double[] aug_img_dy = image.ImageSupplement(image_dy, radius*2+1, radius*2+1);

        Descriptor[] descriptors = new Descriptor[interestPoints.size()];
        for (int k = 0; k < interestPoints.size(); k++)
        {
            descriptors[k] = new Descriptor(barCharCount * basketCount, interestPoints.get(k));
            double[] peaks = getPointOrientation(aug_img_dx, aug_img_dy, interestPoints.get(k), sigma, radius, aug_width);    // Ориентация точки

            for (double phiRotate: peaks)
            {

                InterestPoint.Point p = interestPoints.get(k);
                p.Phi = phiRotate;
                interestPoints.set(k, p);

                for (int i = 1; i < dimension; i++)
                {
                    for (int j = 1; j < dimension; j++)
                    {
                        int x = interestPoints.get(k).x + radius;
                        int y = interestPoints.get(k).y + radius;
                        double gradient_X = aug_img_dx[(j - radius + y)*aug_width + i - radius + x];
                        double gradient_Y = aug_img_dy[(j - radius + y)*aug_width + i - radius + x];

                        // получаем значение(домноженное на Гаусса) и угол
                        double value = getGradientValue(gradient_X, gradient_Y) /* * KernelCreator.getGaussValue(i, j, sigma)*/;
                        double phi = getGradientDirection(gradient_X, gradient_Y) + 2 * Math.PI - phiRotate;
                        phi = (phi % 2 * Math.PI);  // Shift

                        // получаем индекс корзины в которую входит phi и смежную с ней
                        int firstBasketIndex = (int)Math.floor(phi / sector);
                        int secondBasketIndex = (int)(Math.floor((phi - halfSector) / sector) + basketCount) % basketCount;

                        // вычисляем центр
                        double mainBasketPhi = firstBasketIndex * sector + halfSector;

                        // распределяем L(value)
                        double mainBasketValue = (1 - (Math.abs(phi - mainBasketPhi) / sector)) * value;
                        double sideBasketValue = value - mainBasketValue;

                        // вычисляем индекс куда записывать значения
                        int i_Rotate = (int)Math.round((i - radius) * Math.cos(phiRotate) + (j - radius) * Math.sin(phiRotate));
                        int j_Rotate = (int)Math.round(-(i - radius) * Math.sin(phiRotate) + (j - radius) * Math.cos(phiRotate));

                        // отбрасываем
                        if (i_Rotate < -radius || j_Rotate < -radius || i_Rotate >= radius || j_Rotate >= radius)
                        {
                            continue;
                        }

                        int tmp_i = (i_Rotate + radius) / barCharStep;
                        int tmp_j = (j_Rotate + radius) / barCharStep;

                        int indexMain = (tmp_i * barCharCountInLine + tmp_j) * basketCount + firstBasketIndex;
                        int indexSide = (tmp_i * barCharCountInLine + tmp_j) * basketCount + secondBasketIndex;

                        // записываем значения
                        descriptors[k].data[indexMain] += mainBasketValue;
                        descriptors[k].data[indexSide] += sideBasketValue;
                    }
                }
                descriptors[k].normalize();
                descriptors[k].clampData(0, 0.2);
                descriptors[k].normalize();
                descriptors[k].interPoint = interestPoints.get(k);
            }
        }
        return descriptors;
    }

    // Поиск похожих дескрипторов
    public static ArrayList<Vector> findSimilar(Descriptor[] d1, Descriptor[] d2, double treshhold) {
        ArrayList<Vector> similar = new ArrayList<Vector>();
        for (int i = 0; i < d1.length; i++)
        {
            int indexSimilar = -1;
            double prevDistance = Double.MAX_VALUE;       // Предыдущий
            double minDistance = Double.MAX_VALUE;        // Минимальный
            for (int j = 0; j < d2.length; j++)
            {
                double dist = getDistance(d1[i], d2[j]);
                if (dist < minDistance)
                {
                    indexSimilar = j;
                    prevDistance = minDistance;
                    minDistance = dist;
                }
            }

            if (minDistance / prevDistance > treshhold)
            {
                continue;      // отбрасываем
            }
            else
            {
                similar.add(new Vector(d1[i], d2[indexSimilar]));
            }
        }
        return similar;
    }

    private static double getDistance(Descriptor d1, Descriptor d2) {
        double result = 0;
        for (int i = 0; i < d1.data.length; i++)
        {
            double tmp = d1.data[i] - d2.data[i];
            result += tmp * tmp;
        }
        return Math.sqrt(result);
    }
}
