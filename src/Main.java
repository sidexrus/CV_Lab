import ui.Gui;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferInt;
import java.io.File;

public class Main {

    public static void main(String[] args) throws Exception {

        BufferedImage bufimg = ImageIO.read(new File("materials//photo.png"));
        //BufferedImage bufimg = ImageIO.read(new File("materials//sob.jpg"));
        Image img = new Image(bufimg);
        bufimg = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        //double[] kernel = {1,2,1,2,4,2,1,2,1};
        //bufimg.getRaster().setSamples(0,0, img.width, img.height,0, img.Convolution(kernel, 3 , 3));
        //bufimg.getRaster().setSamples(0,0, img.width, img.height,0, img.Sobel());
        bufimg.getRaster().setSamples(0,0, img.getWidth(), img.getHeight(),0, img.GaussianBlur(5,3));
        //bufimg.getRaster().setSamples(0,0, img.width, img.height,0, img.SeparableFilter(new double[]{1,2,1}, new double[]{1,2,1}));
        new Gui().show(img.getImage());
        new Gui().show(bufimg);
    }
}
