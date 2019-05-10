public class Descriptor {
    public InterestPoint.Point interPoint;    // Интересная точка - центр
    public double[] data; // N - Количество корзин * L кол-во гистограмм

    public Descriptor() { }

    public Descriptor(int size, InterestPoint.Point interPoint){
        data = new double[size];
        this.interPoint = interPoint;
    }

    public void normalize(){
        double length = 0;
        for (double a: data)
        length += a * a;

        length = Math.sqrt(length);

        for (int i = 0; i < data.length; i++)
            data[i] /= length;
    }

    public int getSize(){
        return data.length;
    }

    public double getAt(int index) { return data[index]; }

    public InterestPoint.Point getInterPoint() { return interPoint; }

    public static double Clamp(double min, double max, double value)
    {
        if (value < min)
            return min;
        if (value > max)
            return max;
        return value;
    }

    public void clampData(double min, double max){
        for (int i = 0; i < data.length; i++)
            data[i] = Clamp(min, max, data[i]);
    }

    public void setPointXY(int x, int y)
    {
        this.interPoint.x = x;
        this.interPoint.y = y;
    }
}


