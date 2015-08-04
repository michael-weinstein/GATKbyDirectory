import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
class callVariants extends QScript {
	def script() {
		val hc = new HaplotypeCaller
		hc.reference_sequence = new File ("hg42.fa")
		hc.standard_min_confidence_threshold_for_emitting = 10
		hc.standard_min_confidence_threshold_for_calling = 30
		hc.input_file :+= new File ("bogusbams/bogus2.bam")
		hc.out = new File ("bogusbams/bogus2.bam.vcf")
		hc.scatterCount = 20
		hc.memoryLimit = 2
		add(hc)
	}
}
