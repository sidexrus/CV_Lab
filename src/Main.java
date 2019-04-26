import ui.Gui;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferInt;
import java.io.File;

public class Main {

    public static void main(String[] args) throws Exception {

        BufferedImage source = ImageIO.read(new File("materials//photo.png"));
        Image src = new Image(ImageIO.read(new File("materials//photo.png")));
        Image img = new Image(ImageIO.read(new File("materials//photo.png")));
        Image img3 = new Image(ImageIO.read(new File("materials//photo.png")));


        //lab1
        //img.DerivativeX(true);
        //BufferedImage bufimg = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        //bufimg.getRaster().setSamples(0,0, img.getWidth(), img.getHeight(),0, img.Sobel());
        //new Gui().show(img.getImage());
        //new Gui().show(src.getImage());

        //lab2
        //ImagePyramid pyr = new ImagePyramid(src);
        //pyr.CreatePyramid(5,2,0.5, 1);
        //pyr.Save_Pyramid();

        //lab3
        InterestPoint ip = new InterestPoint(img3);
        //new Gui().show(ip.HarrisMap(2));
        new Gui().show(ip.Harris(0.03, 2, 150));
        //new Gui().show(ip.Moravek(0.03, 2, 150));

    }
}
