-- MySQL dump 10.13  Distrib 5.1.51, for pc-linux-gnu (i686)
--
-- Host: 127.0.0.1    Database: cmiData
-- ------------------------------------------------------
-- Server version	5.1.51-debug-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES latin1 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `Study`
--

DROP TABLE IF EXISTS `Study`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Study` (
  `ID` int(8) NOT NULL AUTO_INCREMENT,
  `Name` char(35) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=4 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Study`
--
-- ORDER BY:  `ID`

INSERT INTO `Study` VALUES (1,'Bone Mets');
INSERT INTO `Study` VALUES (2,'CT-Lung');
INSERT INTO `Study` VALUES (3,'Glioma');
INSERT INTO `Study` VALUES (4,'Sarcoma');

--
-- Table structure for table `Species`
--

DROP TABLE IF EXISTS `Species`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Species` (
  `ID` int(1) NOT NULL AUTO_INCREMENT,
  `Name` char(11) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=3 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Species`
--
-- ORDER BY:  `ID`

INSERT INTO `Species` VALUES (1,'Human');
INSERT INTO `Species` VALUES (2,'Rat');
INSERT INTO `Species` VALUES (3,'Mouse');

--
-- Table structure for table `Anatomy`
--

DROP TABLE IF EXISTS `Anatomy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Anatomy` (
  `ID` int(8) NOT NULL AUTO_INCREMENT,
  `Name` char(20) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=3 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Anatomy`
--
-- ORDER BY:  `ID`

INSERT INTO `Anatomy` VALUES (1,'Brain');
INSERT INTO `Anatomy` VALUES (2,'Bone');
INSERT INTO `Anatomy` VALUES (3,'Pancreas');

--
-- Table structure for table `Model`
--

DROP TABLE IF EXISTS `Model`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Model` (
  `ID` int(8) NOT NULL AUTO_INCREMENT,
  `Name` char(20) NOT NULL DEFAULT '',
  `Description` char(100) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Model`
--
-- ORDER BY:  `ID`

INSERT INTO `Model` VALUES (1,'1833','Breast cancer bone metastasis');
INSERT INTO `Model` VALUES (2,'PC3','Prostate cancer bone metastasis');
INSERT INTO `Model` VALUES (3,'9L','Rat gliosarcoma');
INSERT INTO `Model` VALUES (4,'Holland','Breast cancer bone metastasis');
INSERT INTO `Model` VALUES (5,'Sleeping Beauty','Breast cancer bone metastasis');

--
-- Table structure for table `Modality`
--

DROP TABLE IF EXISTS `Modality`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Modality` (
  `ID` int(8) NOT NULL AUTO_INCREMENT,
  `Name` char(20) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=6 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Modality`
--
-- ORDER BY:  `ID`

INSERT INTO `Modality` VALUES (1,'MRI');
INSERT INTO `Modality` VALUES (2,'CT');
INSERT INTO `Modality` VALUES (3,'PET');
INSERT INTO `Modality` VALUES (4,'SPECT');
INSERT INTO `Modality` VALUES (5,'BLI');
INSERT INTO `Modality` VALUES (6,'Histo');

--
-- Table structure for table `Measurement`
--

DROP TABLE IF EXISTS `Measurement`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Measurement` (
  `ID` int(8) NOT NULL AUTO_INCREMENT,
  `Name` char(20) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=2 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Measurement`
--
-- ORDER BY:  `ID`

INSERT INTO `Measurement` VALUES (1,'Volume');
INSERT INTO `Measurement` VALUES (2,'Mean Value');

--
-- Table structure for table `Subject`
--

DROP TABLE IF EXISTS `Subject`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Subject` (
  `ID` int(8) NOT NULL AUTO_INCREMENT,
  `Name` char(20) NOT NULL DEFAULT '',
  `Group` char(10) NOT NULL DEFAULT '',
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=2 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `Subject`
--
-- ORDER BY:  `ID`

INSERT INTO `Subject` VALUES (1,'Volume');
INSERT INTO `Subject` VALUES (2,'Mean Value');

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2010-09-30 11:01:37
